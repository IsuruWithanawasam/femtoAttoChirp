# High Harmonic Generation Chirp Calculator

A Python implementation for calculating attochirp and harmonic chirp rates in High Harmonic Generation (HHG) using the Strong Field Approximation (SFA) and saddle point analysis.

## Overview

This code calculates two types of chirp in HHG spectra:
- **Attochirp (GDD)**: Second-order spectral phase derivative (d²Φ/dΩ²)
- **Harmonic Chirp Rate (b)**: Intensity-dependent phase derivative (dΦ/dI)

The calculations are performed using high-precision arithmetic via `mpmath` to ensure numerical accuracy in solving complex saddle point equations.

## Features

- High-precision calculations (30 decimal places by default)
- Analytical electric field and vector potential functions
- Saddle point equation solver for quantum path analysis
- Quasiclassical action integral computation
- Automatic unit conversion (atomic units ↔ conventional units)
- Short trajectory path analysis

## Requirements

```bash
pip install mpmath
```

## Physical Parameters

### Laser Parameters
- **Wavelength**: 800 nm
- **Pulse Duration**: 8 fs (FWHM)
- **Peak Intensity**: 3×10¹⁴ W/cm²
- **CEP**: 0 (Carrier-Envelope Phase)

### Target Parameters
- **Gas**: Argon
- **Ionization Potential**: 15.76 eV

### Harmonic Analysis
- **Central Harmonic**: 21st order
- **Adjacent Harmonics**: 19th and 23rd (for attochirp calculation)

## Code Structure

### 1. Physical Constants and Unit Conversions
```python
# Conversion factors
FS_TO_AU = 41.341        # 1 fs = 41.341 a.u. of time
EV_TO_AU = 1.0 / 27.2114 # 1 eV = 1/27.2114 a.u.
WCM2_TO_AU = 1.0 / 3.51e16 # Intensity conversion
```

### 2. Core Functions

#### Electric Field
```python
E(t, E0_val)
```
Calculates the time-dependent electric field using a cos² envelope.

#### Vector Potential
```python
A(t, E0_val)
```
Analytical vector potential (negative integral of E-field).

#### Stationary Momentum
```python
p_s(ts, ts_prime, E0_val)
```
Calculates the stationary momentum for a given quantum path.

#### Saddle Point Equations
```python
equations_to_solve(ts, ts_prime, omega_val, E0_val)
```
System of equations F₁=0, F₂=0 that define the saddle points.

#### Quasiclassical Action
```python
get_action_S(ts, ts_prime, ps_val, E0_val)
```
Computes the action integral along the quantum path.

### 3. Main Solver
```python
find_path_phase(harmonic_order_q, intensity_au, initial_guess)
```
Solves saddle point equations and calculates the phase Φ = ωt_s + S.

## Usage

1. **Set precision** (optional):
```python
mp.dps = 30  # Decimal places
```

2. **Configure parameters**:
```python
LAMBDA_NM = mp.mpf('800.0')      # Wavelength (nm)
PULSE_T_FS = mp.mpf('8.0')        # Pulse duration (fs)
INTENSITY_WCM2 = mp.mpf('3e14')   # Peak intensity (W/cm²)
CENTRAL_HARMONIC_ORDER = 21       # Harmonic to analyze
```

3. **Run the script**:
```bash
python hhg_chirp_calculator.py
```

## Output

The code produces the following results:

### Attochirp (GDD)
- Value in atomic units²
- Value in fs²
- Physical interpretation: spectral group delay dispersion

### Harmonic Chirp Rate (b)
- Value in atomic units⁻¹
- Value in fs⁻²
- Physical interpretation: rate of phase change with intensity

### Comparison
- Sign comparison (positive/negative)
- Order of magnitude analysis
- Unit conversions for practical interpretation

## Example Output

```
--- CHIRP RESULTS (Short Path, H21) ---

--- Attochirp (GDD = d2Phi/dOmega^2) ---
  Attochirp = (value) a.u.^2
  Attochirp = (value) fs^2

--- Harmonic Chirp Rate (b) ---
  Harmonic Chirp Rate = (value) a.u.^-1
  Harmonic Chirp Rate = (value) fs^-2

--- Comparison ---
  Sign: Attochirp is positive, Harmonic Chirp is negative.
  Order of Magnitude (a.u.):
    Attochirp (GDD): ~10^X
    Harmonic Chirp Rate (b): ~10^Y
```

## Physics Background

### Strong Field Approximation (SFA)
The code implements the three-step model of HHG:
1. **Ionization**: Electron tunnels through the Coulomb barrier
2. **Propagation**: Electron accelerates in the laser field
3. **Recombination**: Electron recombines, emitting a high-energy photon

### Saddle Point Method
The quantum paths are found by solving saddle point equations that determine:
- **t_s'**: Complex ionization time
- **t_s**: Complex recombination time

### Chirp Types

**Attochirp**: Arises from the intrinsic quantum path dynamics
- Formula: d²Φ/dΩ² ≈ (Φ(q+2) - 2Φ(q) + Φ(q-2)) / (2ω₀)²

**Harmonic Chirp**: Arises from intensity variation during the pulse
- Formula: b ≈ 8ln(2) × I₀/τ² × dΦ/dI

## Mathematical Details

### Phase Calculation
The total phase is:
```
Φ = ω × t_s + S
```
where S is the quasiclassical action:
```
S = ∫[t_s', t_s] [½(p_s - A(t))² + I_p] dt
```

### Jacobian Matrix
For calculating dΦ/dI, the code constructs and inverts the Jacobian:
```
J = [[J₁₁, J₁₂],
     [J₂₁, J₂₂]]
```

## Customization

### Changing the Target Gas
Modify the ionization potential:
```python
IP_EV = mp.mpf('15.76')  # Argon
# IP_EV = mp.mpf('24.59')  # Helium
# IP_EV = mp.mpf('13.60')  # Xenon
```

### Adjusting Laser Parameters
```python
LAMBDA_NM = mp.mpf('800.0')        # Wavelength
PULSE_T_FS = mp.mpf('8.0')          # Pulse duration
INTENSITY_WCM2 = mp.mpf('3e14')    # Peak intensity
CEP = 0                             # Carrier-envelope phase
```

### Selecting Different Harmonics
```python
CENTRAL_HARMONIC_ORDER = 21  # Change to desired harmonic order
```

## Troubleshooting

### Solver Fails
If the saddle point solver fails:
1. Adjust initial guesses (`SHORT_path_guess`)
2. Increase precision (`mp.dps`)
3. Check if parameters are in physical regime

### Numerical Instabilities
- Ensure adequate precision for complex saddle points
- Verify that intensity is not too extreme
- Check that harmonic order is in the physically accessible range

## Author

**Isuru Withanawasam**  
Created: October 20, 2025

## References

This implementation is based on:
- Strong Field Approximation theory
- Lewenstein model for HHG
- Saddle point analysis methods
- Quantum path interference in ultrafast physics

## License

MIT license here

## Citation

If you use this code in your research, please cite appropriately.
