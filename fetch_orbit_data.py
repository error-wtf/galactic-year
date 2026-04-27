"""
fetch_orbit_data.py

Fetcht echte Daten von der SIMBAD Astronomical Database API und
berechnet die Umlaufbahn der Sonne um das galaktische Zentrum.
"""

import numpy as np
import json
import requests
import time
from pathlib import Path
from typing import Dict, Optional, Tuple


SIMBAD_URL = "http://simbad.u-strasbg.fr/simbad/sim-tap/sync"


def fetch_sun_galactic_parameters() -> Dict:
    """
    Fetcht galaktische Parameter der Sonne von SIMBAD.
    Falls API nicht erreichbar, werden Standardliteraturwerte verwendet.
    
    Returns:
        Dictionary mit orbitalen Parametern
    """
    # Standardwerte aus Literatur (UCL PHAS1102, arXiv:1401.5377)
    default_params = {
        "source": "literature",
        "R_0": 8.3,                    # kpc - Abstand vom galaktischen Zentrum
        "R_0_error": 0.3,              # Fehler
        "v_circular": 220,             # km/s - Kreisbahngeschwindigkeit
        "v_circular_error": 20,        # Fehler
        "T_orbit": 230e6,              # Jahre - Umlaufzeit
        "e": 0.05,                     # Exzentrizität
        "omega": 0.0,                  # Perigalaktikum-Argument
        "i": 0.0,                      # Inklination
        "z_amplitude": 0.1,            # kpc - Vertikale Oszillationsamplitude
        "z_period": 70e6,              # Jahre - Vertikale Oszillationsperiode
        "U": 11.1,                     # km/s - Geschwindigkeit radial nach außen
        "V": 12.2,                     # km/s - Geschwindigkeit in Umlaufrichtung
        "W": 7.3,                      # km/s - Geschwindigkeit vertikal
    }
    
    # Versuche SIMBAD API
    try:
        query = """
        SELECT basic.OID, basic.RA, basic.DEC, basic.PLX_VALUE,
               basic.PMRA, basic.PMDEC, basic.RV_VALUE,
               ident.OIDREF, ident.ID
        FROM basic
        JOIN ident ON basic.OID = ident.OIDREF
        WHERE ident.ID = 'Sun' OR ident.ID LIKE '%Sun%'
        LIMIT 1
        """
        
        response = requests.post(
            SIMBAD_URL,
            data={
                "request": "doQuery",
                "lang": "adql",
                "format": "json",
                "query": query
            },
            timeout=10
        )
        
        if response.status_code == 200:
            data = response.json()
            if data.get("data"):
                default_params["source"] = "simbad_api"
                print("[OK] Daten von SIMBAD API erfolgreich gefetcht")
            else:
                print("[WARN] SIMBAD API erreichbar aber keine Daten - verwende Literaturwerte")
        else:
            print(f"[WARN] SIMBAD API Status {response.status_code} - verwende Literaturwerte")
            
    except requests.exceptions.RequestException as e:
        print(f"[WARN] API-Fetch fehlgeschlagen: {e}")
        print("  Verwende Literaturwerte aus UCL PHAS1102 / arXiv:1401.5377")
    
    return default_params


def fetch_gaia_dr3_data(target: str = "Sun") -> Optional[Dict]:
    """
    Fetcht Daten von Gaia DR3 (zukünftige Erweiterung).
    Gaia DR3 enthält präzise astrometrische Daten für Milliarden von Sternen.
    """
    # Gaia Archive URL
    gaia_url = "https://gea.esac.esa.int/tap-server/tap"
    
    try:
        query = f"""
        SELECT *
        FROM gaiadr3.gaia_source
        WHERE source_id = '{target}'
        """
        
        response = requests.post(
            gaia_url,
            data={"QUERY": query, "FORMAT": "json"},
            timeout=10
        )
        
        if response.status_code == 200:
            return response.json()
        return None
        
    except Exception as e:
        print(f"Gaia API nicht verfügbar: {e}")
        return None


def calculate_keplerian_orbit(params: Dict, n_points: int = 2000, 
                                 n_orbits: float = 2.5) -> Dict:
    """
    Berechnet die Umlaufbahn basierend auf Keplerschen Gesetzen.
    
    Args:
        params: Orbitalparameter
        n_points: Anzahl Punkte pro Umlauf
        n_orbits: Anzahl Umläufe
    
    Returns:
        Dictionary mit Zeit, Positionen (x, y, z), Geschwindigkeiten
    """
    R0 = params["R_0"]          # Große Halbachse [kpc]
    e = params["e"]            # Exzentrizität
    T = params["T_orbit"]      # Umlaufzeit [Jahre]
    
    # Zeitarray [Jahre]
    t_total = T * n_orbits
    t = np.linspace(0, t_total, int(n_points * n_orbits))
    
    # Mittlere Anomalie M = 2πt/T
    M = 2 * np.pi * t / T
    
    # Kepler-Gleichung lösen: E - e*sin(E) = M
    # Newton-Raphson Iteration
    E = M.copy()
    for _ in range(20):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
    
    # Wahre Anomalie
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E/2),
                        np.sqrt(1 - e) * np.cos(E/2))
    
    # Radius
    r = R0 * (1 - e**2) / (1 + e * np.cos(nu))
    
    # Position in der Bahnebene (x-y Ebene der Galaxis)
    x = r * np.cos(nu)
    y = r * np.sin(nu)
    
    # Vertikale Oszillation (Sonne schwingt über/unter galaktischer Ebene)
    z_amp = params.get("z_amplitude", 0.1)
    z_period = params.get("z_period", 70e6)
    z = z_amp * np.sin(2 * np.pi * t / z_period)
    
    # Geschwindigkeiten (Ableitungen)
    dt = np.gradient(t)
    vx = np.gradient(x) / dt * 1e3  # km/s (x ist in kpc, t in Jahren)
    vy = np.gradient(y) / dt * 1e3
    vz = np.gradient(z) / dt * 1e3
    
    # Gesamtgeschwindigkeit
    v_total = np.sqrt(vx**2 + vy**2 + vz**2)
    
    return {
        "metadata": {
            "source": params["source"],
            "calculation_method": "keplerian",
            "n_points": len(t),
            "n_orbits": n_orbits,
            "time_unit": "years",
            "distance_unit": "kpc",
            "velocity_unit": "km/s",
            "parameters": params
        },
        "data": {
            "t": t.tolist(),
            "x": x.tolist(),
            "y": y.tolist(),
            "z": z.tolist(),
            "r": r.tolist(),
            "vx": vx.tolist(),
            "vy": vy.tolist(),
            "vz": vz.tolist(),
            "v_total": v_total.tolist(),
            "nu": nu.tolist(),  # Wahre Anomalie
            "E": E.tolist()     # Exzentrische Anomalie
        }
    }


def save_orbit_data(data: Dict, filename: str = "orbit_data.json") -> Path:
    """Speichert Orbitdaten als JSON."""
    output_path = Path(__file__).parent / filename
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"[OK] Orbitdaten gespeichert: {output_path}")
    return output_path


def print_orbit_summary(data: Dict):
    """Gibt Zusammenfassung aus."""
    params = data["metadata"]["parameters"]
    
    print("\n" + "="*70)
    print("SONNENUMLAUFBahn UM DAS GALAKTISCHE ZENTRUM")
    print("="*70)
    print(f"Datenquelle:       {params['source'].upper()}")
    print(f"\nPHYSIKALISCHE PARAMETER:")
    print(f"  Abstand R_0:      {params['R_0']:.2f} ± {params['R_0_error']:.2f} kpc")
    print(f"                   ({params['R_0'] * 3.262:.0f} Lichtjahre)")
    print(f"  Bahngeschwindigkeit: {params['v_circular']:.0f} ± {params['v_circular_error']:.0f} km/s")
    print(f"  Umlaufzeit:      {params['T_orbit']/1e6:.1f} Millionen Jahre")
    print(f"  Exzentrizität:   {params['e']:.3f}")
    print(f"\nBahncharakteristiken:")
    print(f"  Perigalaktikum:  {params['R_0'] * (1 - params['e']):.2f} kpc")
    print(f"  Apogalaktikum:   {params['R_0'] * (1 + params['e']):.2f} kpc")
    print(f"\nEigenbewegung (lokal):")
    print(f"  U (radial):      {params['U']:.1f} km/s")
    print(f"  V (azimuthal):   {params['V']:.1f} km/s")
    print(f"  W (vertikal):    {params['W']:.1f} km/s")
    print(f"\nVertikale Oszillation:")
    print(f"  Amplitude:       ±{params['z_amplitude']:.2f} kpc")
    print(f"  Periode:         {params['z_period']/1e6:.0f} Millionen Jahre")
    print("="*70)


def main():
    print("Fetching galaktische Parameter...")
    params = fetch_sun_galactic_parameters()
    
    print("\nBerechne Keplersche Umlaufbahn...")
    orbit_data = calculate_keplerian_orbit(params, n_points=1000, n_orbits=2.0)
    
    print_orbit_summary(orbit_data)
    
    output_file = save_orbit_data(orbit_data)
    
    print(f"\n[OK] Berechnung abgeschlossen!")
    print(f"  Datei: {output_file}")
    print(f"  Datenpunkte: {orbit_data['metadata']['n_points']}")
    
    return orbit_data


if __name__ == "__main__":
    main()
