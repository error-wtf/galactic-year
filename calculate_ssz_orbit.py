"""
calculate_ssz_orbit.py

Berechnet die Umlaufbahn der Sonne mit Segmented Spacetime (SSZ) Framework.

SSZ-Modifikationen:
- Zeitdilatation D(r) = 1/(1 + Xi(r))
- Weak Field: Xi(r) = r_s/(2r) für r >> r_s
- Strong Field: Xi(r) = 1 - exp(-phi·r/r_s) für r ~ r_s
- phi = (1 + √5)/2 ≈ 1.618 (Goldener Schnitt)

Die SSZ-Zeit beeinflusst die effektive Umlaufgeschwindigkeit und Periode.
"""

import numpy as np
import json
from pathlib import Path
from typing import Dict, Tuple
import sys


# Goldener Schnitt
PHI = (1 + np.sqrt(5)) / 2  # ≈ 1.618033988749895

# Physikalische Konstanten
G = 6.67430e-11           # m³ kg⁻¹ s⁻²
M_sun = 1.98847e30        # kg
M_galactic = 4.3e6 * M_sun  # Masse Sgr A* (4.3 Millionen Sonnenmassen)
c = 299792458             # m/s
AU = 1.496e11             # m
kpc_to_m = 3.086e19       # 1 kpc in Metern
year_to_s = 365.25 * 24 * 3600  # Jahr in Sekunden


def schwarzschild_radius(M: float) -> float:
    """Berechnet Schwarzschild-Radius: r_s = 2GM/c²"""
    return 2 * G * M / c**2


def xi_weak_field(r: np.ndarray, r_s: float) -> np.ndarray:
    """
    SSZ Weak Field Approximation: Xi(r) = r_s/(2r)
    Gültig für r/r_s >> 100
    """
    return r_s / (2 * r)


def xi_strong_field(r: np.ndarray, r_s: float, phi: float = PHI) -> np.ndarray:
    """
    SSZ Strong Field: Xi(r) = 1 - exp(-phi·r/r_s)
    Gültig für r/r_s < 100
    """
    return 1 - np.exp(-phi * r / r_s)


def xi_unified(r: np.ndarray, r_s: float, r_transition: float = 100.0, 
                phi: float = PHI) -> np.ndarray:
    """
    Unified SSZ Xi(r) mit glattem Übergang zwischen weak und strong field.
    
    Args:
        r: Abstand in Einheiten von r_s (r/r_s)
        r_s: Schwarzschild-Radius
        r_transition: Übergangsbereich (default 100)
        phi: Goldener Schnitt
    
    Returns:
        Xi(r) als ndarray
    """
    # Berechne beide Formen
    xi_weak = xi_weak_field(r, r_s)
    xi_strong = xi_strong_field(r, r_s, phi)
    
    # Glatter Übergang mit Sigmoid-Funktion
    # Für r >> r_transition: ξ_weak
    # Für r << r_transition: ξ_strong
    x = np.log10(r / r_s)  # Logarithmischer Abstand
    x0 = np.log10(r_transition)
    
    # Sigmoid-Übergang
    sigmoid = 1 / (1 + np.exp(-2 * (x - x0)))
    
    # Kombiniere
    xi = xi_weak * sigmoid + xi_strong * (1 - sigmoid)
    
    return xi


def ssz_time_dilation(xi: np.ndarray) -> np.ndarray:
    """
    SSZ Zeitdilatation: D(r) = 1/(1 + Xi(r))
    
    Dies ist der entscheidende Unterschied zu GR:
    - GR: Zeitdilatation geht gegen 0 bei r → r_s
    - SSZ: Zeitdilatation erreicht D(r_s) ≈ 0.556 (finit!)
    """
    return 1.0 / (1.0 + xi)


def ssz_orbital_velocity(r: float, M: float, xi: float) -> float:
    """
    Berechnet Orbitalgeschwindigkeit mit SSZ-Korrektur.
    
    Klassisch: v = √(GM/r)
    SSZ: v_ssz = v_klassisch · D(r) = √(GM/r) / (1 + Xi(r))
    
    Die Zeitdilatation verringert die effektive Geschwindigkeit.
    """
    v_classical = np.sqrt(G * M / r)
    D = ssz_time_dilation(np.array([xi]))[0]
    return v_classical * D


def calculate_ssz_orbit(params: Dict, n_points: int = 2000, 
                         n_orbits: float = 2.0) -> Dict:
    """
    Berechnet Umlaufbahn mit SSZ-Korrekturen.
    
    Args:
        params: Orbitalparameter
        n_points: Anzahl Punkte
        n_orbits: Anzahl Umläufe
    
    Returns:
        Dictionary mit klassischen und SSZ-Orbitdaten
    """
    # Parameter
    R0_kpc = params["R_0"]  # kpc
    R0_m = R0_kpc * kpc_to_m
    e = params["e"]
    T_classical = params["T_orbit"]  # Jahre
    
    # Galaktische Masse (Sgr A*)
    M_bh = 4.3e6 * M_sun  # 4.3 Millionen Sonnenmassen
    r_s = schwarzschild_radius(M_bh)
    
    print(f"\n  Schwarzschild-Radius Sgr A*: {r_s/1e3:.2f} km = {r_s:.2e} m")
    print(f"  Sonnenabstand: {R0_kpc:.2f} kpc = {R0_m:.2e} m")
    print(f"  r/r_s = {R0_m/r_s:.2e}")
    
    # Zeitarray
    t_total = T_classical * n_orbits
    t = np.linspace(0, t_total, int(n_points * n_orbits))
    
    # Kepler-Bahn (klassisch)
    M = 2 * np.pi * t / T_classical
    E = M.copy()
    for _ in range(20):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
    
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2),
                        np.sqrt(1 - e) * np.cos(E/2))
    
    # Radius in Metern
    r_m = R0_m * (1 - e**2) / (1 + e * np.cos(nu))
    
    # SSZ Xi(r) für jeden Punkt
    xi_values = xi_unified(r_m, r_s)
    
    # SSZ Zeitdilatation
    D_values = ssz_time_dilation(xi_values)
    
    # Klassische Geschwindigkeit
    v_classical = np.sqrt(G * M_bh / r_m)  # m/s
    
    # SSZ-korrigierte Geschwindigkeit
    v_ssz = v_classical * D_values
    
    # SSZ-korrigierte Periode (Zeitdilatation verlängert effektive Periode)
    # T_ssz = ∫ dt / D(r(t)) - effektive Periode ist länger
    T_ssz_factor = np.mean(1.0 / D_values)
    T_ssz = T_classical * T_ssz_factor
    
    # Positionen in kpc
    x_kpc = (r_m * np.cos(nu)) / kpc_to_m
    y_kpc = (r_m * np.sin(nu)) / kpc_to_m
    
    # Vertikale Oszillation
    z_amp = params.get("z_amplitude", 0.1) * kpc_to_m
    z_period = params.get("z_period", 70e6) * year_to_s
    t_s = t * year_to_s
    z_m = z_amp * np.sin(2 * np.pi * t_s / z_period)
    z_kpc = z_m / kpc_to_m
    
    return {
        "metadata": {
            "calculation_method": "ssz",
            "n_points": len(t),
            "n_orbits": n_orbits,
            "time_unit": "years",
            "distance_unit": "kpc",
            "velocity_unit": "m/s",
            "parameters": params,
            "ssz_parameters": {
                "phi": float(PHI),
                "M_bh_solar_masses": 4.3e6,
                "r_s_meters": float(r_s),
                "r_s_km": float(r_s/1000)
            },
            "comparison": {
                "T_classical_years": float(T_classical),
                "T_ssz_years": float(T_ssz),
                "T_factor_ssz": float(T_ssz_factor),
                "v_mean_classical_km_s": float(np.mean(v_classical)/1000),
                "v_mean_ssz_km_s": float(np.mean(v_ssz)/1000),
                "xi_at_R0": float(xi_values[len(xi_values)//2]),
                "D_at_R0": float(D_values[len(D_values)//2])
            }
        },
        "data": {
            "t": t.tolist(),
            "x": x_kpc.tolist(),
            "y": y_kpc.tolist(),
            "z": z_kpc.tolist(),
            "r": (r_m / kpc_to_m).tolist(),
            "xi": xi_values.tolist(),
            "D": D_values.tolist(),
            "v_classical": v_classical.tolist(),
            "v_ssz": v_ssz.tolist(),
            "nu": nu.tolist()
        }
    }


def print_ssz_summary(data: Dict):
    """Gibt SSZ-Zusammenfassung aus."""
    meta = data["metadata"]
    comp = meta["comparison"]
    ssz_params = meta["ssz_parameters"]
    
    print("\n" + "="*70)
    print("SSZ (SEGMENTED SPACETIME) BERECHNUNG")
    print("="*70)
    print(f"Goldener Schnitt phi = {ssz_params['phi']:.10f}")
    print(f"\nSCHWARZSCHILD-RADIUS Sgr A*:")
    print(f"  r_s = {ssz_params['r_s_km']:.2f} km = {ssz_params['r_s_meters']:.2e} m")
    print(f"\nSSZ ZEITDILATATION bei R_0:")
    print(f"  Xi(R_0) = {comp['xi_at_R0']:.8f}")
    print(f"  D(R_0) = 1/(1+Xi) = {comp['D_at_R0']:.8f}")
    print(f"\nGESCHWINDIGKEITSVERGLEICH:")
    print(f"  Klassisch:  {comp['v_mean_classical_km_s']:.1f} km/s")
    print(f"  SSZ:        {comp['v_mean_ssz_km_s']:.1f} km/s")
    print(f"  Reduktion:  {(1 - comp['v_mean_ssz_km_s']/comp['v_mean_classical_km_s'])*100:.4f}%")
    print(f"\nUMLAUFZEITVERGLEICH:")
    print(f"  Klassisch:  {comp['T_classical_years']/1e6:.1f} Millionen Jahre")
    print(f"  SSZ:        {comp['T_ssz_years']/1e6:.1f} Millionen Jahre")
    print(f"  Faktor:     {comp['T_factor_ssz']:.6f}x länger")
    print("="*70)


def save_ssz_data(data: Dict, filename: str = "ssz_orbit_data.json") -> Path:
    """Speichert SSZ-Daten."""
    output_path = Path(__file__).parent / filename
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"[OK] SSZ-Daten gespeichert: {output_path}")
    return output_path


def main():
    print("SSZ Orbitberechnung")
    print("=" * 50)
    
    # Lade Parameter (oder verwende Defaults)
    params = {
        "source": "ssz_calculation",
        "R_0": 8.3,          # kpc
        "R_0_error": 0.3,
        "v_circular": 220,   # km/s
        "v_circular_error": 20,
        "T_orbit": 230e6,    # Jahre
        "e": 0.05,
        "z_amplitude": 0.1,  # kpc
        "z_period": 70e6     # Jahre
    }
    
    print("Berechne SSZ-korrigierte Umlaufbahn...")
    ssz_data = calculate_ssz_orbit(params, n_points=1000, n_orbits=2.0)
    
    print_ssz_summary(ssz_data)
    
    output_file = save_ssz_data(ssz_data)
    
    print(f"\n[OK] SSZ-Berechnung abgeschlossen!")
    print(f"  Datei: {output_file}")
    
    return ssz_data


if __name__ == "__main__":
    main()
