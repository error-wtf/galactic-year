import requests, numpy as np

# Test SIMBAD URLs
for url in [
    "https://simbad.cds.unistra.fr/simbad/tap/sync",
    "https://simbad.u-strasbg.fr/simbad/tap/sync",
    "http://simbad.u-strasbg.fr/simbad/tap/sync",
    "https://simbad.cds.unistra.fr/simbad/tap/",
]:
    adql = "SELECT TOP 1 main_id, pmra, pmdec FROM basic WHERE main_id = 'NAME Sgr A'"
    try:
        r = requests.post(url, data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"json","QUERY":adql}, timeout=10)
        print(f"{url}  ->  {r.status_code}")
        if r.status_code == 200:
            d = r.json()
            cols = [c["name"] for c in d["metadata"]]
            rows = [dict(zip(cols, row)) for row in d["data"]]
            print(f"   DATA: {rows}")
            break
    except Exception as e:
        print(f"{url}  ->  ERROR: {e}")

# Check Parallaxen: 0.306 mas fuer HMSFRs
# Reid+2014 Table 1 - typische Parallaxen sind 0.1-0.5 mas fuer ~2-10 kpc Distanz
plx_median_mas = 0.306
dist_kpc = 1.0 / (plx_median_mas / 1000.0)   # FALSCH: /1000 macht aus mas arcsec -> kpc
dist_kpc2 = 1000.0 / plx_median_mas           # RICHTIG: 1kpc = 1/1arcsec, 1arcsec=1000mas
print(f"\nParallaxe 0.306 mas:")
print(f"  1.0 / (plx/1000) = {dist_kpc:.1f} kpc  (FALSCH - umrechnung fehlerhaft)")
print(f"  1000.0 / plx     = {dist_kpc2:.1f} kpc  (RICHTIG)")
print(f"\nTypische HMSFR bei 3.27 kpc: plx = {1000/3.27:.3f} mas  (korrekt fuer W51 etc.)")
