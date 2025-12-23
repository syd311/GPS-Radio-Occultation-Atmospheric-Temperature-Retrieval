# GPS-Radio-Occultation-Atmospheric-Temperature-Retrieval

## Project Overview
This project implements an end-to-end Remote Sensing Retrieval Pipeline for GPS Radio Occultation (GPS-RO) data. It simulates the physics of the **COSMIC-2** satellite mission to retrieve atmospheric temperature profiles from space.

The pipeline performs two key functions:
1.  **Forward Modeling:** Simulates satellite "Bending Angle" telemetry based on a Standard Atmosphere model (including realistic sensor noise).
2.  **Inverse Retrieval:** Solves the Inverse Abel Transform and integrates the Hydrostatic Equation to recover the temperature profile from the noisy telemetry.

---

## Results (Validation)
The algorithm successfully retrieves the temperature profile (Blue), matching the "True" input atmosphere (Red) within **<1% error** in the core region (10km–40km). It correctly identifies the **Tropopause** at ~11km.

![Validation Plot](/atmospheric_retrieval/results/result.png)

*(Note: This plot is generated automatically by running the pipeline)*

---

This project solves the classic **Inverse Problem** of Radio Occultation.

### 1. The Abel Inversion
The core challenge is converting the satellite's raw measurement (Bending Angle, $\alpha$) into Atmospheric Refractivity ($n$). This requires solving a Volterra integral equation of the second kind numerically:

$$
n(a) = \exp\left( \frac{1}{\pi} \int_{a}^{\infty} \frac{\alpha(x)}{\sqrt{x^2 - a^2}} \, dx \right)
$$

*   **Implementation:** Discretized using the Trapezoidal rule with singularity handling at the lower bound ($x=a$).

### 2. Hydrostatic Integration
Once Refractivity ($N$) is known, density is derived. Pressure ($P$) is calculated by integrating the air column from the Top of Atmosphere (TOA) downwards:

$$
P(z) = \int_{z}^{z_{top}} \rho(h) \cdot g(h) \, dh
$$

Temperature ($T$) is then retrieved via the Ideal Gas Law relation for Refractivity:

$$
T(z) = 77.6 \cdot \frac{P(z)}{N(z)}
$$

---
