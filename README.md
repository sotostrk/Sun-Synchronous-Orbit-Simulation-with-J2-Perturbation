# Sun-Synchronous Orbit Simulation with J2 Perturbation and Solar Direction Tracking

This project simulates a satellite in a Sun-synchronous low Earth orbit (SSO), incorporating the secular effects of the J2 zonal harmonic perturbation due to Earth's oblateness. The orbit is propagated in the Earth-Centered Inertial (ECI) frame using a modified two-body dynamic model. 

The simulation visualizes how the J2 effect causes a steady drift in the Right Ascension of the Ascending Node (RAAN), enabling the orbit plane to maintain a constant angle relative to the Sun. The solar direction is modeled as a rotating unit vector to approximate Earth's revolution around the Sun.
The result is a physically accurate demonstration of Sun-synchronous orbital behavior with 3D visualization and animation.

---

##  Objectives

- Simulate the precessional effect of Earth’s J2 zonal harmonic on orbital RAAN.
- Compute Sun-synchronous inclinations analytically for a given orbital altitude.
- Visualize the orbital trajectory and Sun vector evolution in the ECI frame.
- Provide both static and animated views to support interpretation of RAAN drift.

---

##  Background

Sun-synchronous orbits are designed such that the orbital plane precesses at a rate equal to Earth’s mean motion around the Sun (~0.9856°/day). 
This ensures the satellite crosses any given point on Earth at the same local solar time. The precession is driven primarily by the oblateness of the Earth, modeled here using the J2 perturbation.


---

## Simulation Features

- **Orbital dynamics with J2 perturbation**: Implements standard first-order zonal harmonic terms in the acceleration model.
- **Analytical inclination computation**: Ensures Sun-synchronous RAAN precession rate given the orbital semi-major axis.
- **Sun direction modeling**: Approximates solar motion in ECI using a circular ecliptic orbit.
- **RAAN drift analysis**: Extracts and plots RAAN angle over time to verify precession rate.
- **Static and animated visualization**: Full 3D orbital trail with Earth and a Sun vector rendered via `matplotlib`.

---

###  Orbit Animation (Download to View)

This MP4 shows the evolving Sun-synchronous orbit with J2-induced RAAN drift and a rotating solar vector.

 [Click to download and view the animation](https://github.com/sotostrk/Sun-Synchronous-Orbit-Simulation-with-J2-Perturbation/raw/main/sun_sync_orbit.mp4)

---

##  Output Examples

### 3D Orbit Plot with Sun Vector
Displays the full orbital trajectory and a snapshot of the solar direction vector at mid-simulation.

![3D Orbit Over 60 Days](https://raw.githubusercontent.com/sotostrk/Sun-Synchronous-Orbit-Simulation-with-J2-Perturbation/main/orbit_60days.png)


---

### RAAN Drift Plot
Demonstrates the secular RAAN change due to the J2 perturbation.

![RAAN Drift Over Time](https://raw.githubusercontent.com/sotostrk/Sun-Synchronous-Orbit-Simulation-with-J2-Perturbation/main/raan_drift.png)










