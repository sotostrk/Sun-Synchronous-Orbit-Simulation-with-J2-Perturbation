import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# parameters
mu_earth = 398600.0  # km^3/s^2
earth_radius = 6378.0  # km
Re = earth_radius
J2 = 1.08263e-3         

# converting orbital elements to position and velocity vectors in ECI frame
def coe_to_rv(a, e, i_deg, raan_deg, argp_deg, nu_deg, mu):
    i = np.radians(i_deg)
    raan = np.radians(raan_deg)
    argp = np.radians(argp_deg)
    nu = np.radians(nu_deg)

    # Compute distance and speed in perifocal frame
    p = a * (1 - e**2)
    r_mag = p / (1 + e * np.cos(nu))

    # Position and velocity in perifocal coordinates (PQW)
    r_pqw = np.array([
        r_mag * np.cos(nu),
        r_mag * np.sin(nu),
        0.0
    ])
    v_pqw = np.array([
        -np.sqrt(mu / p) * np.sin(nu),
        np.sqrt(mu / p) * (e + np.cos(nu)),
        0.0
    ])

    # Rotation matrix to convert PQW to ECI
    Rz_raan = np.array([
        [np.cos(raan), -np.sin(raan), 0],
        [np.sin(raan), np.cos(raan), 0],
        [0, 0, 1]
    ])
    Rx_i = np.array([
        [1, 0, 0],
        [0, np.cos(i), -np.sin(i)],
        [0, np.sin(i), np.cos(i)]
    ])
    Rz_argp = np.array([
        [np.cos(argp), -np.sin(argp), 0],
        [np.sin(argp), np.cos(argp), 0],
        [0, 0, 1]
    ])

    # Final rotation matrix
    R = Rz_raan @ Rx_i @ Rz_argp

    # Convert to ECI frame
    r_eci = R @ r_pqw
    v_eci = R @ v_pqw

    return r_eci, v_eci
#

def sun_vector_eci(t_seconds):
    # Earth orbital angular velocity around the Sun
    days_per_year = 365.25
    omega_sun = 2 * np.pi / (days_per_year * 86400)  # rad/s

    theta = omega_sun * t_seconds  # angle from x-axis
    r_sun = np.array([
        np.cos(theta),
        np.sin(theta),
        0.0
    ])
    return r_sun  # Unit vector in ECI frame

#  two-body gravitational dynamics
def two_body_dynamics_j2(t, y):
    rx, ry, rz, vx, vy, vz = y
    r = np.array([rx, ry, rz])
    norm_r = np.linalg.norm(r)

    # Basic two-body gravitational acceleration
    acc_gravity = -mu_earth * r / norm_r**3

    # J2 perturbation
    z2 = rz**2
    r2 = norm_r**2
    factor = 1.5 * J2 * mu_earth * Re**2 / norm_r**5

    ax_j2 = factor * rx * (1 - 5 * z2 / r2)
    ay_j2 = factor * ry * (1 - 5 * z2 / r2)
    az_j2 = factor * rz * (3 - 5 * z2 / r2)
    acc_j2 = np.array([ax_j2, ay_j2, az_j2])

    acc_total = acc_gravity + acc_j2
    return [vx, vy, vz, acc_total[0], acc_total[1], acc_total[2]]

#

sim_days = 60
#  orbit in 3D
def plot_trajectory(r, times):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    #  orbit trail
    ax.plot(r[:, 0], r[:, 1], r[:, 2], color='m', linewidth=0.6, label='Orbit')
    ax.plot([r[0, 0]], [r[0, 1]], [r[0, 2]], 'ok', label='Start')

    #  Earth
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)
    ax.plot_surface(x, y, z, cmap='Blues', alpha=0.7)

    # Sun vector 
    mid_index = len(times) // 2
    sun_vec = sun_vector_eci(times[mid_index])

    ax.quiver(0, 0, 0,
              sun_vec[0], sun_vec[1], sun_vec[2],
              color='orange', length=earth_radius * 2,
              normalize=True, label='Sun Vector')

    #view limits
    max_val = np.max(np.abs(r)) * 1.1
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.set_title(f'Orbit Over {sim_days} Days with J2 Drift')
    ax.legend()
    plt.show()



def sun_sync_inclination(a_km):
    J2 = 1.08263e-3
    Re = 6378.0  # km
    mu = 398600.0  # km^3/s^2

    # Target RAAN rate for Sun-synchronous orbit: 360°/365.25 days = 0.9856 deg/day
    dOmega_dt_target = -0.9856 * (np.pi / 180) / 86400  # rad/s

    factor = -2 * dOmega_dt_target / (3 * J2 * (Re / a_km)**2 * np.sqrt(mu / a_km**3))
    if abs(factor) > 1.0:
        raise ValueError("Invalid orbit altitude — no real solution for inclination.")
    return np.degrees(np.arccos(factor))


#
e = 0
perigee_altitude = 1050  # km above Earth surface
r_p = earth_radius + perigee_altitude
a = r_p / (1 - e)  # Set semi-major axis to keep perigee realistic
i = sun_sync_inclination(a)  #  Sun-synchronous inclination
raan = 0.0               # RAAN (deg)
argp = 0.0                # Argument of perigee (deg) initial value
nu = 0.0                  # True anomaly at start (deg)

# orbital elements to state vector
r0, v0 = coe_to_rv(a, e, i, raan, argp, nu, mu_earth)
y0 = np.hstack((r0, v0))

# Time span for one full orbit
orbital_period = 2 * np.pi * np.sqrt(a**3 / mu_earth)
# Simulate for 3 months
t_final = sim_days * 86400  # seconds
t_span = (0, t_final)
t_eval = np.linspace(t_span[0], t_span[1], int(sim_days * 96))  
# 96 points per day around 1 point every 15 minutes

# Integrate 
sol = solve_ivp(two_body_dynamics_j2, t_span, y0, t_eval=t_eval, method='RK45', rtol=1e-9, atol=1e-9)

##
velocity = sol.y[3:6, :].T  # Satellite velocity over time
times = sol.t

sun_angles = []
##
r_eci = sol.y[:3, :].T
v_eci = sol.y[3:6, :].T
raan_list = []

for i in range(len(times)):
    r = r_eci[i]
    v = v_eci[i]

    h = np.cross(r, v)
    n = np.cross([0, 0, 1], h)

    if np.linalg.norm(n) < 1e-6:
        raan = 0.0
    else:
        n_unit = n / np.linalg.norm(n)
        raan = np.arccos(np.clip(n_unit[0], -1.0, 1.0))
        if n_unit[1] < 0:
            raan = 2 * np.pi - raan
    raan_list.append(np.degrees(raan))
    
    
    #plot actual raan drift 
plt.figure()
plt.plot(times / 86400, raan_list)
plt.xlabel('Time (days)')
plt.ylabel('RAAN (deg)')
plt.title('RAAN Drift Over Time Due to J2')
plt.grid()
plt.pause(0.001)

for i in range(len(times)):
    v = velocity[i]
    sun_vec = sun_vector_eci(times[i])
    
    # Normalize vectors
    v_unit = v / np.linalg.norm(v)
    s_unit = sun_vec / np.linalg.norm(sun_vec)
    
    # Angle between them (in degrees)
    dot = np.clip(np.dot(v_unit, s_unit), -1.0, 1.0)
    angle_deg = np.degrees(np.arccos(dot))
    sun_angles.append(angle_deg)




# Extract and plot trajectory
r_traj = sol.y[:3, :].T
plot_trajectory(r_traj, times)





# === LIVE ANIMATION FUNCTION ===
def animate_orbit_gif(r_traj_local, times_local):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Plot Earth
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)
    ax.plot_surface(x, y, z, cmap='Blues', alpha=0.5)

    # Set plot limits
    max_val = np.max(np.linalg.norm(r_traj_local, axis=1)) * 1.1
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    ax.set_zlabel("Z (km)")
    ax.set_title("Sun-Synchronous Orbit with J2 Drift")

    trail, = ax.plot([], [], [], 'm', lw=1.5, label='Orbit')
    sun_arrow = [None]

    def update(frame):
        r_now = r_traj_local[:frame + 1]
        trail.set_data(r_now[:, 0], r_now[:, 1])
        trail.set_3d_properties(r_now[:, 2])

        if sun_arrow[0] is not None:
            sun_arrow[0].remove()

        sun_vec = sun_vector_eci(times_local[frame])
        sun_arrow[0] = ax.quiver(0, 0, 0,
                                 sun_vec[0], sun_vec[1], sun_vec[2],
                                 color='orange', length=earth_radius * 2, normalize=True)

        return trail, sun_arrow[0]

    ani = animation.FuncAnimation(fig, update, frames=len(r_traj_local), interval=20, blit=False)

    
    ani.save("sun_sync_orbit.mp4", writer='ffmpeg', fps=30, bitrate=1800)
    

plt.legend()
plt.show()

if __name__ == "__main__":
    animate_orbit_gif(r_traj, times)
