# Cavity Analysis 🏮

Analysis of the transmon-cavity system used in our lab.

## 👉 Half vs Full

Comparison of a half Tesla cavity to a full Tesla cavity. Checking for quality
factors, life-times and dissipation using both a volume integral method
(modeling the dirt as a 3D object) and surface integral method (integrting the
surface of the cavity to analyze the dirt).

## 👉 Coupling

Analysis of the coupling between the transom and the readout strip-line.

**Note:** For this to work you need pyEPR version from after 19/05/2020 (see [this commit](https://github.com/zlatko-minev/pyEPR/commit/29ee909e1fa2c0ddd879afac4ca90123098b8baf)).

## 👉 Chip Test

GDS design of chip used to test the cavity with a coil instead of a real transmon. "Fake" (linear) transmon-readout chip on a saphire substrate.

## 👉 Dephasing

Simulations and analytical solutions of the dephasing problem, used to imporve the dephasing time of the cavity in the transmon-cavity system.
