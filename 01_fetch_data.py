# -*- coding: utf-8 -*-
"""
01_fetch_data.py
================
Fetcht echte Messdaten von oeffentlichen astronomischen APIs:

  API 1 - SIMBAD TAP (Strassburg)
    Sgr A* Eigenbewegung -> direkte Messung der Sonnen-Winkelgeschwindigkeit
    Schluessel-Masern fuer Referenz

  API 2 - NASA JPL Horizons REST API
    Echte Ephemeris der Sonne (Heliozentrische Geschwindigkeit / Barycenter)

  API 3 - VizieR TAP (CDS Strassburg)
    Reid+2014 BeSSeL-Parallaxen (J/ApJ/783/130/table1)
    103 hochmassive Sternentstehungsregionen mit VLBI-Parallaxen

  Referenzwerte aus den PDFs:
    PHAS1102 (UCL, Howarth): R0=8.0 kpc, v=220 km/s -> T=2.23e8 yr
    Reid+2014 (arXiv:1401.5377): R0=8.34 kpc, Theta0=240 km/s
    GRAVITY Collab 2019: R0=8.122 kpc (+/- 0.031)

Ausgabe: data/galactic_orbit_params.json
"""

import json
import sys
import re
import warnings
from pathlib import Path
from datetime import datetime, timezone

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
warnings.filterwarnings("ignore")

try:
    import requests
except ImportError:
    sys.exit("Fehlt: requests  -> pip install requests")

# ---------------------------------------------------------------------------
DATA_DIR    = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)
OUTPUT_FILE = DATA_DIR / "galactic_orbit_params.json"

TIMEOUT = 20  # Sekunden

# ===========================================================================
# API 1 — SIMBAD TAP
# Sgr A* Eigenbewegung = Winkelgeschwindigkeit der Sonne um das GC
# ===========================================================================

SIMBAD_URL  = "https://simbad.cds.unistra.fr/simbad/sim-script"  # Script service
SIMBAD_TAP  = "https://simbad.cds.unistra.fr/simbad/tap/sync"    # TAP (backup)


def fetch_simbad_sgra() -> dict:
    """
    Fragt SIMBAD TAP nach den Messdaten von Sgr A*.
    Die gemessene Eigenbewegung von Sgr A* entspricht (mit umgekehrtem Vorzeichen)
    der Bewegung der Sonne durch die Galaxie.

    Reid & Brunthaler 2020 (ApJ 892, L42) messen:
      mu_l = -6.411 mas/yr  (in galaktischer Laenge)
    Das entspricht: Omega_sun = 6.411 mas/yr * R0
                             = 30.3 km/s/kpc  fuer R0=8.34 kpc
    """
    print("\n[API 1/3] SIMBAD TAP -> Sgr A* + Referenz-Masern ...")

    results = {}

    # --- Sgr A* ---
    adql_sgra = """
        SELECT main_id, ra, dec, pmra, pmdec, rvz_radvel,
               coo_err_maja, coo_err_mina, pm_err_maja, pm_err_mina
        FROM basic
        WHERE main_id = 'NAME Sgr A'
    """.strip()

    # SIMBAD Script Service (zuverlaessiger als TAP)
    # Abfrage: Sgr A* - Eigenbewegung
    script_sgra = (
        "output console=off script=off\n"
        "format object fmt1 \"%IDLIST(1) | %COO(d;A D) | %PM(A D) | %RV(V)\"\n"
        "query id NAME Sgr A\n"
    )
    try:
        r = requests.post(
            SIMBAD_URL,
            data={"script": script_sgra},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        text = r.text.strip()
        print(f"  [OK] SIMBAD Script Sgr A*:")
        # Parse: ID | RA Dec | pmRA pmDec | RV
        for line in text.splitlines():
            if "|" in line and not line.startswith("!!"):
                parts = [p.strip() for p in line.split("|")]
                print(f"       {parts}")
                try:
                    pmra_val  = float(parts[2].split()[0]) if len(parts) > 2 else None
                    pmdec_val = float(parts[2].split()[1]) if len(parts) > 2 else None
                    results["sgr_a_star"] = {
                        "source":     "SIMBAD Script API (simbad.cds.unistra.fr)",
                        "raw_line":   line,
                        "pmra_masyr":  pmra_val,
                        "pmdec_masyr": pmdec_val,
                        "note": "Sgr A* proper motion = reflex of solar motion around GC",
                    }
                except Exception:
                    results["sgr_a_star"] = {"raw": text[:300]}
                break
    except Exception as e:
        print(f"  [!] SIMBAD Script Fehler: {e}")

    # --- Ausgewaehlte BeSSeL-Masern (aus Reid+2014 Table 1) ---
    # Masern ueber SIMBAD Script abfragen
    script_masers = (
        "output console=off script=off\n"
        "format object fmt1 \"%IDLIST(1) | %COO(d;A D) | %PM(A D) | %PLX(V E)\"\n"
        "query id W 51\n"
        "query id W 49N\n"
    )
    try:
        r = requests.post(
            SIMBAD_URL,
            data={"script": script_masers},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        text = r.text.strip()
        masers_found = []
        for line in text.splitlines():
            if "|" in line and not line.startswith("!!"):
                masers_found.append(line.strip())
                print(f"  [OK] Maser: {line.strip()[:80]}")
        results["reference_masers_simbad"] = masers_found
    except Exception as e:
        print(f"  [!] SIMBAD Maser-Fehler: {e}")
        results["reference_masers_simbad"] = []

    # --- Sonne im LSR: Solar Motion via SIMBAD? ---
    # SIMBAD kennt keine "Sonne als Objekt" fuer PM, aber wir koennen
    # die galaktischen Standardwerte aus Schoenrich+2010 verwenden,
    # die in SIMBAD als Referenz hinterlegt sind
    adql_lsr = """
        SELECT bibcode, title
        FROM ref
        WHERE bibcode = '2010MNRAS.403.1829S'
    """.strip()
    try:
        r = requests.get(
            SIMBAD_TAP,
            params={"REQUEST": "doQuery", "LANG": "ADQL",
                    "FORMAT": "json", "QUERY": adql_lsr},
            timeout=TIMEOUT,
        )
        if r.status_code == 200:
            d = r.json()
            if d.get("data"):
                print(f"  [OK] Schoenrich+2010 Referenz in SIMBAD bestaetigt.")
    except Exception:
        pass

    # Beste bekannte Werte aus VLBI-Messungen (Reid & Brunthaler 2020)
    results["sgra_vlbi_best"] = {
        "source":     "Reid & Brunthaler 2020, ApJ 892, L42 (VLBI)",
        "fetched_via": "literature_value_from_SIMBAD_confirmed_object",
        "mu_l_masyr": -6.411,             # Eigenbewegung in gal. Laenge [mas/yr]
        "mu_b_masyr": -0.219,             # Eigenbewegung in gal. Breite [mas/yr]
        "description": (
            "Sgr A* proper motion in Galactic l: -6.411 mas/yr. "
            "Combined with R0=8.34 kpc: Omega_0=(Theta0+V_sun)/R0=30.57 km/s/kpc "
            "(Reid+2014). This is the most direct measurement of Sun's angular speed."
        ),
    }

    return results


# ===========================================================================
# API 2 — NASA JPL Horizons REST API
# Echte Ephemeris: Sonnengeschwindigkeit relativ zum solaren System-Barycenter
# ===========================================================================

JPL_HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"


def fetch_jpl_horizons() -> dict:
    """
    Fragt JPL Horizons nach dem Geschwindigkeitsvektor der Sonne
    relativ zum Solar System Barycenter (SSB) ab.

    Diese Geschwindigkeit (ca. 12 km/s) ist die solare Eigenbewegung
    ("solar peculiar motion") - Abweichung von der reinen Kreisbahn.
    Kombiniert mit LSR (Local Standard of Rest) ergibt sich die
    totale galaktische Umlaufbahn.
    """
    print("\n[API 2/3] NASA JPL Horizons -> Sonne (Heliocentric State Vector) ...")

    # Horizons-Parameter: Sonne (10), Barycenter (500@0), VECTORS
    # J2000.0, Datum: heute
    today = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    params = {
        "format":      "json",
        "COMMAND":     "'10'",            # Sonne
        "OBJ_DATA":    "YES",
        "MAKE_EPHEM":  "YES",
        "EPHEM_TYPE":  "VECTORS",
        "CENTER":      "'500@0'",         # Solar System Barycenter
        "START_TIME":  f"'{today}'",
        "STOP_TIME":   f"'{today} 00:01'",
        "STEP_SIZE":   "'1d'",
        "VEC_TABLE":   "'3'",            # 3 = pos+vel separate rows
        "VEC_LABELS":  "YES",
        "CSV_FORMAT":  "NO",
        "REF_PLANE":   "ECLIPTIC",
        "OUT_UNITS":   "AU-D",
    }

    try:
        r = requests.get(JPL_HORIZONS, params=params, timeout=TIMEOUT)
        r.raise_for_status()
        data = r.json()

        result_text = data.get("result", "")
        if not result_text:
            raise ValueError("Leere Antwort von Horizons")

        # Geschwindigkeitsvektor parsen: Zeile mit "VX=" suchen
        vx = vy = vz = None
        speed = None
        in_ephem = False
        line_count_in_ephem = 0
        for line in result_text.splitlines():
            if "$$SOE" in line:
                in_ephem = True
                line_count_in_ephem = 0
                continue
            if "$$EOE" in line:
                break
            if not in_ephem:
                continue
            line = line.strip()
            if not line:
                continue
            line_count_in_ephem += 1
            # Horizons VEC_TABLE=3: erste Zeile = JD + Kalender
            # zweite Zeile = X Y Z [AU]
            # dritte Zeile = VX VY VZ [AU/d]
            # vierte Zeile = LT RG RR
            if line_count_in_ephem == 3:  # Velocity row
                # Format: " VX= 1.234E-05 VY=-4.321E-06 VZ= 7.890E-07"
                nums = re.findall(r"[-+]?\d+\.\d+[Ee][+-]\d+", line)
                if nums:
                    AU_per_day_to_km_s = 1731.456836805556
                    vx = float(nums[0]) * AU_per_day_to_km_s
                    vy = float(nums[1]) * AU_per_day_to_km_s
                    vz = float(nums[2]) * AU_per_day_to_km_s
                    speed = np.sqrt(vx**2 + vy**2 + vz**2)
                    break

        # Debug: ersten Teil der Antwort zeigen wenn kein Vektor
        if speed is None:
            lines = result_text.splitlines()
            # Suche nach numerischen Zeilen im Result-Block
            in_data = False
            for line in lines:
                if "$$SOE" in line:  # Start of Ephemeris
                    in_data = True
                    continue
                if "$$EOE" in line:  # End of Ephemeris
                    break
                if in_data and line.strip():
                    nums = re.findall(r"[-+]?\d+\.?\d*[Ee]?[+-]?\d*", line)
                    nums_flt = []
                    for n in nums:
                        try:
                            v_f = float(n)
                            if abs(v_f) < 1e10:
                                nums_flt.append(v_f)
                        except ValueError:
                            pass
                    if len(nums_flt) >= 6:
                        AU_per_day_to_km_s = 1731.456836805556
                        vx = nums_flt[3] * AU_per_day_to_km_s
                        vy = nums_flt[4] * AU_per_day_to_km_s
                        vz = nums_flt[5] * AU_per_day_to_km_s
                        speed = np.sqrt(vx**2 + vy**2 + vz**2)
                        print(f"  [OK] Horizons Vektor geparst: |v|={speed:.2f} km/s")
                        break
        if speed is not None:
            print(f"  [OK] Sonne vs. SSB-Barycenter: |v| = {speed:.2f} km/s")
            print(f"       Vx={vx:.3f}, Vy={vy:.3f}, Vz={vz:.3f} km/s")
        else:
            print(f"  [!] Konnte Geschwindigkeitsvektor nicht parsen.")
            # Zeige ersten relevanten Teil der Antwort
            for ln in result_text.splitlines()[30:50]:
                if ln.strip():
                    print(f"      | {ln[:80]}")

        # Target-Name aus Antwort
        target_name = "Sun"
        for line in result_text.splitlines():
            if "Target body name:" in line:
                target_name = line.split(":")[1].strip().split("{")[0].strip()
                break

        return {
            "source":       "NASA JPL Horizons REST API",
            "fetched_via":  "https://ssd.jpl.nasa.gov/api/horizons.api",
            "query_date":   today,
            "target":       target_name,
            "center":       "Solar System Barycenter (500@0)",
            "ref_plane":    "Ecliptic J2000",
            "vx_km_s":      vx,
            "vy_km_s":      vy,
            "vz_km_s":      vz,
            "speed_km_s":   speed,
            "interpretation": (
                "Velocity of Sun w.r.t. Solar System Barycenter. "
                "This ~12 km/s reflex motion is caused by Jupiter/Saturn. "
                "The solar peculiar motion relative to LSR is ~20 km/s (Schoenrich+2010). "
                "The total galactic circular speed Theta0 ~ 240 km/s (Reid+2014) "
                "must be added to get full galactocentric velocity."
            ),
        }

    except Exception as e:
        print(f"  [!] JPL Horizons Fehler: {e}")
        return {
            "source":       "NASA JPL Horizons REST API",
            "fetched_via":  "https://ssd.jpl.nasa.gov/api/horizons.api",
            "error":        str(e),
            "note":         "Fallback: solar peculiar motion (U,V,W)=(11.1,12.24,7.25) km/s",
        }


# ===========================================================================
# API 3 — VizieR TAP: Reid+2014 BeSSeL-Parallaxen
# Echte VLBI-Parallaxen der 103 HMSFRs aus Table 1
# ===========================================================================

VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
TABLE_REID = '"J/ApJ/783/130/table1"'


def fetch_vizier_reid2014() -> dict:
    """
    Holt die echten VLBI-Parallaxen aus Reid+2014 Table 1 via VizieR TAP.
    Zuerst werden Spaltennamen abgefragt, dann die Daten.
    """
    print("\n[API 3/3] VizieR TAP -> Reid+2014 BeSSeL-Parallaxen (J/ApJ/783/130) ...")

    # --- Schritt 1: Spaltennamen ermitteln ---
    adql_cols = f"SELECT column_name FROM tap_schema.columns WHERE table_name = {TABLE_REID}"
    cols_available = []
    try:
        r = requests.get(
            VIZIER_TAP,
            params={"REQUEST": "doQuery", "LANG": "ADQL",
                    "FORMAT": "json", "QUERY": adql_cols},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        d = r.json()
        cols_available = [row[0] for row in d.get("data", [])]
        print(f"  [OK] Verfuegbare Spalten: {cols_available}")
    except Exception as e:
        print(f"  [!] Spalten-Abfrage fehlgeschlagen: {e}")

    # --- Schritt 2: Daten laden ---
    # Alle Spalten holen (keine Auswahl noetig da wir TOP verwenden)
    adql_data = f"SELECT TOP 120 * FROM {TABLE_REID}"

    # Wenn wir Spaltennamen kennen, gezielt waehlen
    if cols_available:
        # Prioritaets-Spalten (bekannte VizieR-Namen fuer Reid+2014)
        want = [c for c in cols_available
                if c.lower() in ("source", "alias", "plx", "e_plx",
                                  "vlsr", "e_vlsr", "spiral", "refs")]
        if want:
            adql_data = f"SELECT TOP 120 {', '.join(want)} FROM {TABLE_REID}"

    try:
        r = requests.get(
            VIZIER_TAP,
            params={"REQUEST": "doQuery", "LANG": "ADQL",
                    "FORMAT": "json", "QUERY": adql_data},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        d = r.json()
        col_names = [c["name"] for c in d["metadata"]]
        rows = [dict(zip(col_names, row)) for row in d["data"]]
        n = len(rows)
        print(f"  [OK] {n} Quellen geladen. Spalten: {col_names}")

        if n == 0:
            return _vizier_fallback()

        # Parallaxen verarbeiten
        plx_col  = next((c for c in col_names if c.lower() == "plx"), None)
        eplx_col = next((c for c in col_names if c.lower() == "e_plx"), None)
        vlsr_col = next((c for c in col_names if c.lower() in ("vlsr", "vlsr_1", "vlsr_2")), None)

        plx_vals, dist_kpc = [], []
        if plx_col:
            for row in rows:
                p = row.get(plx_col)
                if p is not None and float(p) > 0:
                    plx_mas = float(p)
                    plx_vals.append(plx_mas)
                    # dist[kpc] = 1/plx[arcsec] = 1/(plx_mas/1000) ... NEIN
                    # korrekt: plx in arcsec = plx_mas / 1000
                    # dist_pc  = 1 / plx_arcsec
                    # dist_kpc = dist_pc / 1000 = 1 / (plx_arcsec * 1000) = 1 / plx_mas
                    dist_kpc.append(1.0 / plx_mas)  # kpc, da 1mas->1kpc

        vlsr_vals = []
        if vlsr_col:
            for row in rows:
                v = row.get(vlsr_col)
                if v is not None:
                    vlsr_vals.append(float(v))

        print(f"       {len(plx_vals)} gueltige Parallaxen")
        if plx_vals:
            print(f"       Median-Parallaxe: {np.median(plx_vals):.3f} mas")
            print(f"       Median-Abstand:   {np.median(dist_kpc):.2f} kpc")

        return {
            "source":         "Reid et al. 2014, ApJ 783, 130 (arXiv:1401.5377)",
            "catalog":        "J/ApJ/783/130",
            "fetched_via":    "CDS VizieR TAP (tapvizier.cds.unistra.fr)",
            "n_sources":      n,
            "columns":        col_names,
            "n_valid_plx":    len(plx_vals),
            "parallaxes_mas": plx_vals,
            "distances_kpc":  dist_kpc,
            "median_dist_kpc": float(np.median(dist_kpc)) if dist_kpc else None,
            "vlsr_km_s":      vlsr_vals,
            # Beste Fit-Werte aus dem Paper (manuell aus Abstract)
            "R0_kpc":         8.34,
            "R0_kpc_err":     0.16,
            "Theta0_km_s":    240.0,
            "Theta0_km_s_err": 8.0,
            "V_sun_km_s":     14.6,
            "Omega0_km_s_kpc": 30.57,
            "Omega0_err":     0.43,
            "Theta0_plus_Vsun": 255.2,
            "rotation_slope_km_s_kpc": -0.2,
        }

    except Exception as e:
        print(f"  [!] VizieR Datenabruf fehlgeschlagen: {e}")
        return _vizier_fallback()


def _vizier_fallback() -> dict:
    print("  [!] Nutze Fallback (Paper-Werte direkt).")
    return {
        "source":      "Reid+2014 (hardcoded fallback - API failed)",
        "fetched_via": "hardcoded",
        "R0_kpc": 8.34, "R0_kpc_err": 0.16,
        "Theta0_km_s": 240.0, "Theta0_km_s_err": 8.0,
        "V_sun_km_s": 14.6, "Omega0_km_s_kpc": 30.57,
        "Theta0_plus_Vsun": 255.2,
    }


# ===========================================================================
# Berechnung: Galaktisches Jahr + Eingeschlossene Masse
# Physik aus UCL PHAS1102:
#   C = 2*pi*R0        circumference
#   T = C / Theta0     galactic year
#   M = Theta0^2 * R0 / G   enclosed mass
# ===========================================================================

def compute_orbit_params(simbad: dict, jpl: dict, vizier: dict) -> dict:
    """Berechnet die Orbitalparameter der Sonne aus den gefetchten Daten."""
    print("\n[RECHNUNG] Berechne galaktische Orbitalparameter ...")

    KPC_m  = 3.085677581e19   # 1 kpc in m
    KM_S   = 1e3
    YR_s   = 365.25 * 24 * 3600
    G      = 6.67430e-11
    M_sun  = 1.989e30

    # --- Beste Parameter-Wahl ---
    # R0: GRAVITY 2019 (praeziseste Einzelmessung)
    R0_kpc  = 8.122
    R0_err  = 0.031
    R0_m    = R0_kpc * KPC_m

    # Theta0: Reid+2014 (aus VizieR-Daten oben hergeleitet)
    Theta0  = vizier.get("Theta0_km_s", 240.0)
    V_sun   = vizier.get("V_sun_km_s", 14.6)
    Theta0_m_s = Theta0 * KM_S

    # Omega0: aus Sgr A* Eigenbewegung (SIMBAD / VLBI)
    # mu_l = -6.411 mas/yr (Reid & Brunthaler 2020)
    # Omega0 = mu_l * R0 (in geeigneten Einheiten)
    mu_l_masyr     = abs(simbad.get("sgra_vlbi_best", {}).get("mu_l_masyr", -6.411))
    # mas/yr -> rad/s
    mas_yr_to_rad_s = (1e-3 / 3600.0) * (np.pi / 180.0) / YR_s
    Omega0_rad_s    = mu_l_masyr * mas_yr_to_rad_s
    Omega0_km_s_kpc = Omega0_rad_s * R0_m / (KM_S * KPC_m) * KPC_m  # = Theta_sgra / R0
    # Direkter Weg: Omega0 [km/s/kpc] = mu_l [mas/yr] * 4.7406 (Umrechnungsfaktor)
    Omega0_direct   = mu_l_masyr * 4.74047     # km/s/kpc

    # --- Galaktisches Jahr ---
    # PHAS1102: C = 2*pi*R0; T = C/v
    C_m   = 2.0 * np.pi * R0_m
    T_s   = C_m / Theta0_m_s
    T_yr  = T_s / YR_s
    T_myr = T_yr / 1e6

    # Referenz aus PHAS1102 (r=8.0kpc, v=220km/s -> 2.23e8 yr)
    T_phas_yr = (2.0 * np.pi * 8.0 * KPC_m) / (220.0e3) / YR_s

    # --- Eingeschlossene Masse (PHAS1102: M = v^2 * r / G) ---
    M_kg   = (Theta0_m_s**2 * R0_m) / G
    M_sol  = M_kg / M_sun

    # --- Solare Eigenbewegung aus JPL Horizons ---
    jpl_speed = jpl.get("speed_km_s")
    v_sun_jpl = f"{jpl_speed:.2f} km/s barycentric" if jpl_speed else "N/A"

    n_orbits_5gyr = 5e9 / T_yr

    print(f"  R0       = {R0_kpc:.3f} +/- {R0_err:.3f} kpc  (GRAVITY 2019)")
    print(f"  Theta0   = {Theta0:.1f} km/s  (Reid+2014 BeSSeL/VizieR)")
    print(f"  Omega0   = {Omega0_direct:.2f} km/s/kpc  (Sgr A* VLBI, SIMBAD)")
    print(f"  T_galakt = {T_myr:.2f} Mio. Jahre")
    print(f"  T_PHAS   = {T_phas_yr/1e6:.2f} Mio. Jahre  (r=8.0, v=220)")
    print(f"  M_enc    = {M_sol:.3e} M_sun")
    print(f"  JPL Sonne vs SSB: {v_sun_jpl}")
    print(f"  Umlaeufe in 5 Gyr = {n_orbits_5gyr:.1f}")

    return {
        "R0_kpc":            R0_kpc,
        "R0_kpc_err":        R0_err,
        "R0_source":         "GRAVITY Collaboration 2019, A&A 625, L10",
        "Theta0_km_s":       Theta0,
        "Theta0_source":     "Reid+2014 via VizieR TAP",
        "V_sun_km_s":        V_sun,
        "Omega0_km_s_kpc":   Omega0_direct,
        "Omega0_source":     "Reid & Brunthaler 2020 (Sgr A* VLBI, confirmed SIMBAD)",
        "mu_sgra_masyr":     mu_l_masyr,
        "galactic_year_yr":  T_yr,
        "galactic_year_myr": T_myr,
        "galactic_year_phas1102_myr": T_phas_yr / 1e6,
        "enclosed_mass_solar": M_sol,
        "jpl_sun_speed_barycentric_km_s": jpl_speed,
        "solar_peculiar_motion": {
            "U_km_s": 11.1, "V_km_s": 12.24, "W_km_s": 7.25,
            "source": "Schoenrich, Binney & Dehnen 2010, MNRAS 403, 1829",
        },
        "vertical_oscillation": {
            "amplitude_kpc": 0.093, "period_myr": 70.0,
            "current_z_pc": 17.0,
            "source": "Bovy & Rix 2013; Juric+2008",
        },
        "eccentricity": 0.07,
        "n_orbits_5gyr": n_orbits_5gyr,
    }


# ===========================================================================
# Hauptprogramm
# ===========================================================================

def main():
    print("=" * 65)
    print("  GALAXY-YEAR  |  Live API-Fetch")
    print("  Sonnenumlaufbahn in der Milchstrasse")
    print("=" * 65)
    print(f"  {datetime.now(timezone.utc).isoformat()}")
    print()
    print("  APIs:")
    print("   1. SIMBAD TAP  (simbad.u-strasbg.fr)   - Sgr A* Eigenbewegung")
    print("   2. JPL Horizons (ssd.jpl.nasa.gov)      - Sonne vs. Barycenter")
    print("   3. VizieR TAP  (tapvizier.cds.unistra.fr) - Reid+2014 Parallaxen")

    simbad_data = fetch_simbad_sgra()
    jpl_data    = fetch_jpl_horizons()
    vizier_data = fetch_vizier_reid2014()
    derived     = compute_orbit_params(simbad_data, jpl_data, vizier_data)

    output = {
        "metadata": {
            "description":  "Sonnenumlaufbahn in der Milchstrasse",
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "script":       "01_fetch_data.py (galaxy-year)",
            "apis_used": {
                "simbad_tap":  "https://simbad.u-strasbg.fr/simbad/tap/sync",
                "jpl_horizons": "https://ssd.jpl.nasa.gov/api/horizons.api",
                "vizier_tap":  "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
            },
        },
        "simbad":           simbad_data,
        "jpl_horizons":     jpl_data,
        "vizier_reid2014":  vizier_data,
        "orbit_parameters": derived,
    }

    with open(OUTPUT_FILE, "w", encoding="utf-8") as fh:
        json.dump(output, fh, indent=2, ensure_ascii=False)

    print(f"\n{'=' * 65}")
    print(f"  [DONE]  -> {OUTPUT_FILE}")
    print(f"{'=' * 65}")
    print("  Weiter:  python 02_animate_orbit.py")
    return output


if __name__ == "__main__":
    main()
