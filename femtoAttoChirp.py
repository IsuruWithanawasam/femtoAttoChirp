# -*- coding: utf-8 -*-
"""
Created on Mon Oct 20 23:35:59 2025

@author: Isuru Withanawasam
"""

import mpmath as mp

# --- 1. Set Precision (Decimal Places) ---
mp.dps = 30

print(f"Using mpmath with {mp.dps} decimal places of precision.")

# --- 2. Define Physical Constants and Parameters ---
# --- Laser Parameters ---
LAMBDA_NM = mp.mpf('800.0')      # Wavelength (nm)
PULSE_T_FS = mp.mpf('8.0')        # Pulse duration (fs) - THIS IS FWHM
INTENSITY_WCM2 = mp.mpf('3e14') # Peak intensity (W/cm^2)
CEP = 0                           # Carrier-Envelope Phase (psi)

# --- Target Parameters (Argon) ---
IP_EV = mp.mpf('15.76')           # Ionization potential of Argon (eV)

# --- Harmonic to Solve For (Central Harmonic) ---
CENTRAL_HARMONIC_ORDER = 21

# --- Conversion to Atomic Units (a.u.) ---
FS_TO_AU = mp.mpf('41.341')       # 1 fs = 41.341 a.u. of time
EV_TO_AU = mp.mpf(1.0 / 27.2114)  # 1 eV = 1/27.2114 a.u. (Hartree)
WCM2_TO_AU = mp.mpf(1.0 / 3.51e16) # 1 a.u. intensity = 3.51e16 W/cm^2

# --- Final Parameters in a.u. ---
Ip = IP_EV * EV_TO_AU                             # Ionization potential (a.u.)
w0 = mp.mpf('45.563') / LAMBDA_NM                  # Central frequency (a.u.)
T0 = 2.0 * mp.pi / w0                             # Optical cycle (a.u.)
I_central = INTENSITY_WCM2 * WCM2_TO_AU           # Central Intensity (a.u.)
T_pulse_fwhm_au = PULSE_T_FS * FS_TO_AU           # Pulse duration FWHM (a.u.)
# tau from paper's pulse definition T = 2*tau*arccos(2^(-1/4)) - envelope parameter, not FWHM
tau_envelope = T_pulse_fwhm_au / (2.0 * mp.acos(mp.power(2, -0.25)))

print(f"--- Central Atomic Units ---")
print(f"Ip (a.u.): {Ip}")
print(f"w0 (a.u.): {w0} (T0 = {T0} a.u.)")
print(f"I_central (a.u.): {I_central}")
print(f"Pulse FWHM (a.u.): {T_pulse_fwhm_au}")
print(f"tau_envelope (a.u.): {tau_envelope}")
print("-" * 20)

# --- 3. Define Functions (need to accept E0) ---
w1 = w0
w2 = w0 + 2.0 / tau_envelope
w3 = w0 - 2.0 / tau_envelope

def E(t, E0_val):
    """ Analytical Electric Field E(t) """
    # Use tau_envelope for the pulse shape
    return E0_val * (mp.cos(t/tau_envelope))**2 * mp.cos(w0*t + CEP)

def A(t, E0_val):
    """ Analytical Vector Potential A(t) = -integral(E(t)) """
    term1 = mp.sin(w1 * t + CEP) / w1
    term2 = 0.5 * mp.sin(w2 * t + CEP) / w2
    term3 = 0.5 * mp.sin(w3 * t + CEP) / w3
    return -E0_val / 2.0 * (term1 + term2 + term3)

def S_A(t, E0_val):
    """ Analytical Integral of A(t), S_A(t) = integral(A(t)) """
    term1 = -mp.cos(w1 * t + CEP) / (w1**2)
    term2 = -0.5 * mp.cos(w2 * t + CEP) / (w2**2)
    term3 = -0.5 * mp.cos(w3 * t + CEP) / (w3**2)
    return -E0_val / 2.0 * (term1 + term2 + term3)

def p_s(ts, ts_prime, E0_val):
    """ Calculates the stationary momentum p_s """
    if ts == ts_prime:
        return mp.mpc('1e100', '1e100')
    integral_A = S_A(ts, E0_val) - S_A(ts_prime, E0_val)
    return (1.0 / (ts - ts_prime)) * integral_A

def equations_to_solve(ts, ts_prime, omega_val, E0_val):
    """ System of equations F1=0, F2=0 """
    try:
        ps_val = p_s(ts, ts_prime, E0_val)
        A_ts = A(ts, E0_val)
        A_ts_prime = A(ts_prime, E0_val)
    except (ZeroDivisionError, OverflowError):
        return (mp.mpc('1e100', '1e100'), mp.mpc('1e100', '1e100'))
    F1 = 0.5 * (ps_val - A_ts_prime)**2 + Ip
    F2 = omega_val - 0.5 * (ps_val - A_ts)**2 - Ip
    return (F1, F2)

def get_action_S(ts, ts_prime, ps_val, E0_val):
    """ Calculates the quasiclassical action S """
    integrand = lambda t: 0.5 * (ps_val - A(t, E0_val))**2 + Ip
    S_val = mp.quad(integrand, [ts_prime, ts])
    return S_val

def calculate_phase_and_deriv_I(ts, ts_prime, omega_val, Intensity_val, E0_val):
    """ Calculates Phi = omega*ts + S and dPhi/dI """
    # This function remains mostly the same as before
    ps_val = p_s(ts, ts_prime, E0_val)
    A_ts = A(ts, E0_val)
    A_ts_prime = A(ts_prime, E0_val)
    E_ts = E(ts, E0_val)
    E_ts_prime = E(ts_prime, E0_val)

    S_val = get_action_S(ts, ts_prime, ps_val, E0_val)
    Phi = omega_val * ts + S_val # Using Phi = omega*ts + S
    dS_dI_partial = (S_val - Ip * (ts - ts_prime)) / Intensity_val

    dp_dts = (A_ts - ps_val) / (ts - ts_prime)
    dp_dts_prime = (-A_ts_prime + ps_val) / (ts - ts_prime)

    J_11 = (ps_val - A_ts_prime) * dp_dts
    J_12 = (ps_val - A_ts_prime) * (dp_dts_prime + E_ts_prime)
    J_21 = -(ps_val - A_ts) * (dp_dts + E_ts)
    J_22 = -(ps_val - A_ts) * dp_dts_prime
    J = mp.matrix([[J_11, J_12], [J_21, J_22]])

    try:
        J_inv = J**-1
    except ZeroDivisionError:
        print(f"Warning: Jacobian inversion failed for I={Intensity_val}. Cannot calculate dPhi/dI.")
        dPhi_dI_nan = mp.mpc(mp.nan, mp.nan)
        return Phi, dPhi_dI_nan

    dps_dI = ps_val / (2 * Intensity_val)
    dAts_dI = A_ts / (2 * Intensity_val)
    dAts_prime_dI = A_ts_prime / (2 * Intensity_val)
    dF1_dI = (ps_val - A_ts_prime) * (dps_dI - dAts_prime_dI)
    dF2_dI = -(ps_val - A_ts) * (dps_dI - dAts_dI)
    dF_dI_vec = mp.matrix([dF1_dI, dF2_dI])
    dtimes_dI_vec = - J_inv * dF_dI_vec
    dts_dI = dtimes_dI_vec[0]
    dPhi_dI = 2 * omega_val * dts_dI + dS_dI_partial

    return Phi, dPhi_dI

def calculate_phase_only(ts, ts_prime, omega_val, E0_val):
    """ Calculates only the phase Phi = omega*ts + S """
    ps_val = p_s(ts, ts_prime, E0_val)
    S_val = get_action_S(ts, ts_prime, ps_val, E0_val)
    Phi = omega_val * ts + S_val
    return Phi

# --- 4. Main Solver Function ---
def find_path_phase(harmonic_order_q, intensity_au, initial_guess):
    """
    Solves saddle point equations for a given harmonic and intensity,
    then calculates Phi = omega*ts + S.
    Returns (ts, ts_prime, Phi) or (None, None, None) on failure.
    """
    omega_val = harmonic_order_q * w0
    E0_val = mp.sqrt(intensity_au)
    try:
        ts_sol, ts_prime_sol = mp.findroot(
            lambda ts, ts_prime: equations_to_solve(ts, ts_prime, omega_val, E0_val),
            (initial_guess[0], initial_guess[1])
        )
        Phi_sol = calculate_phase_only(ts_sol, ts_prime_sol, omega_val, E0_val)
        
        # *** FIXED LINE 169 ***
        print(f"  Solver success for H{harmonic_order_q} at I={intensity_au}") 
        
        return ts_sol, ts_prime_sol, Phi_sol
    except Exception as e:
        
        # *** FIXED LINE 172 ***
        print(f"  Solver failed for H{harmonic_order_q} at I={intensity_au}: {e}") 
        
        return None, None, None

# --- 5. Define Initial Guesses ---
# Using the guess that worked for H21 short path
ts_prime_re_guess_short = mp.mpf('25.0')
ts_re_guess_short = mp.mpf('70.0')
SHORT_path_guess = (
    mp.mpc(ts_re_guess_short, mp.mpf('1.5')),
    mp.mpc(ts_prime_re_guess_short, mp.mpf('1.5'))
)

# --- 6. Calculate Phase for q-2, q, q+2 (Short Path Only) ---
q_minus_2 = CENTRAL_HARMONIC_ORDER - 2 # H19
q_central = CENTRAL_HARMONIC_ORDER     # H21
q_plus_2  = CENTRAL_HARMONIC_ORDER + 2 # H23

print("\n--- Calculating Phase for Attochirp (Short Path) ---")
ts_s_m2, ts_p_s_m2, Phi_s_m2 = find_path_phase(q_minus_2, I_central, SHORT_path_guess)
ts_s_c,  ts_p_s_c,  Phi_s_c  = find_path_phase(q_central, I_central, SHORT_path_guess)
ts_s_p2, ts_p_s_p2, Phi_s_p2 = find_path_phase(q_plus_2,  I_central, SHORT_path_guess)

# --- 7. Calculate dPhi/dI for q (Short Path Only) ---
print("\n--- Calculating dPhi/dI for Harmonic Chirp (Short Path, H{}) ---".format(q_central))
try:
    # We already have ts_s_c, ts_p_s_c from the step above
    if ts_s_c is not None:
         # *** THIS LINE IS FIXED: Use the correct omega for q_central ***
         omega_central = q_central * w0 
         _, dPhi_dI_s_c = calculate_phase_and_deriv_I(ts_s_c, ts_p_s_c, omega_central, I_central, mp.sqrt(I_central))
    else:
         dPhi_dI_s_c = None
         print("  Cannot calculate dPhi/dI because central solver failed.")
except Exception as e:
    dPhi_dI_s_c = None
    print(f"  Calculation of dPhi/dI failed: {e}")


# --- 8. Calculate Chirps ---

# Attochirp (d2Phi/dOmega^2)
attochirp = mp.mpc(mp.nan, mp.nan) # Default to NaN
delta_omega = 2 * w0
if Phi_s_m2 is not None and Phi_s_c is not None and Phi_s_p2 is not None:
    attochirp = (Phi_s_p2 - 2*Phi_s_c + Phi_s_m2) / (delta_omega**2)

# Harmonic Chirp Rate (b)
harmonic_chirp_rate = mp.mpc(mp.nan, mp.nan) # Default to NaN
if dPhi_dI_s_c is not None:
    # Using formula from slides: b ~ 8*ln(2) * I0 / tau_FWHM^2 * dPhi/dI
    harmonic_chirp_rate = 8 * mp.log(2) * I_central / (T_pulse_fwhm_au**2) * dPhi_dI_s_c


# --- 9. Display Final Results ---
print("\n" + "="*40)
print(f"      CHIRP RESULTS (Short Path, H{CENTRAL_HARMONIC_ORDER})")
print("="*40 + "\n")

print(f"--- Attochirp (GDD = d2Phi/dOmega^2) ---")
if mp.isnan(attochirp.real):
    print("  Calculation failed (likely due to solver failure for adjacent harmonics).")
else:
    print(f"  Attochirp = {attochirp} a.u.^2")
    # Convert to fs^2 (1 a.u. time = 24.19 as = 0.02419 fs)
    attochirp_fs2 = attochirp * (mp.mpf('0.0241888')**2)
    print(f"  Attochirp = {attochirp_fs2} fs^2")


print(f"\n--- Harmonic Chirp Rate (b) ---")
if mp.isnan(harmonic_chirp_rate.real):
     print("  Calculation failed (likely due to dPhi/dI failure).")
else:
    print(f"  Harmonic Chirp Rate = {harmonic_chirp_rate} a.u.^-1")
    # Convert to fs^-2 (1 a.u. time = 24.19 as = 0.02419 fs)
    harmonic_chirp_rate_fs_inv2 = harmonic_chirp_rate / (FS_TO_AU**2)
    print(f"  Harmonic Chirp Rate = {harmonic_chirp_rate_fs_inv2} fs^-2")


print("\n--- Comparison ---")
if not mp.isnan(attochirp.real) and not mp.isnan(harmonic_chirp_rate.real):
    # Sign
    sign_atto = "positive" if attochirp.real > 0 else "negative"
    sign_harmonic = "positive" if harmonic_chirp_rate.real > 0 else "negative"
    print(f"  Sign: Attochirp is {sign_atto}, Harmonic Chirp is {sign_harmonic}.")
    
    # Order of Magnitude (using fs^-2 for comparison)
    # Convert attochirp (fs^2) to chirp rate (fs^-2) for rough comparison
    # This isn't a direct comparison of the same physical quantity, but relates them
    # A simple estimate: Rate ~ 1 / Duration^2; Duration ~ sqrt(GDD)
    # This is very approximate
    
    # Let's compare the raw atomic unit values first
    log10_atto = mp.log10(mp.fabs(attochirp.real)) if attochirp.real != 0 else -mp.inf
    log10_harmonic = mp.log10(mp.fabs(harmonic_chirp_rate.real)) if harmonic_chirp_rate.real != 0 else -mp.inf
    
    print(f"  Order of Magnitude (a.u.):")
    print(f"    Attochirp (GDD): ~10^{int(mp.floor(log10_atto))}")
    print(f"    Harmonic Chirp Rate (b): ~10^{int(mp.floor(log10_harmonic))}")
    
    # Rough comparison in fs^-2
    print(f"  Order of Magnitude (fs^-2):")
    log10_atto_fs2 = mp.log10(mp.fabs(attochirp_fs2.real)) if attochirp_fs2.real != 0 else -mp.inf
    log10_harmonic_fs_inv2 = mp.log10(mp.fabs(harmonic_chirp_rate_fs_inv2.real)) if harmonic_chirp_rate_fs_inv2.real != 0 else -mp.inf
    print(f"    Attochirp (GDD in fs^2): ~10^{int(mp.floor(log10_atto_fs2))}")
    print(f"    Harmonic Chirp Rate (b in fs^-2): ~10^{int(mp.floor(log10_harmonic_fs_inv2))}")
else:
    print("  Cannot compare magnitudes because one or both calculations failed.")