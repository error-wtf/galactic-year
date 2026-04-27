"""
run_all.py

Master-Script: Fetched Daten, berechnet klassisch + SSZ, erstellt alle Animationen.
"""

import subprocess
import sys
from pathlib import Path


def run_script(script_name: str, description: str) -> bool:
    """Führt ein Python-Script aus."""
    print(f"\n{'='*70}")
    print(f"  {description}")
    print(f"  Script: {script_name}")
    print(f"{'='*70}")
    
    script_path = Path(__file__).parent / script_name
    if not script_path.exists():
        print(f"  [WARN] Script nicht gefunden: {script_path}")
        return False
    
    try:
        result = subprocess.run(
            [sys.executable, str(script_path)],
            cwd=script_path.parent,
            check=True,
            capture_output=False
        )
        print(f"  [OK] {script_name} erfolgreich ausgeführt")
        return True
    except subprocess.CalledProcessError as e:
        print(f"  [ERR] Fehler in {script_name}: {e}")
        return False


def main():
    print("\n" + "="*70)
    print("  GALACTIC YEAR - VOLLSTÄNDIGE PIPELINE")
    print("="*70)
    
    # Erstelle Output-Verzeichnis
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)
    
    scripts = [
        ("fetch_orbit_data.py", "SCHRITT 1: Fetch Orbitdaten (API + Literatur)"),
        ("calculate_ssz_orbit.py", "SCHRITT 2: SSZ-Berechnung (Segmented Spacetime)"),
        ("animate_orbit.py", "SCHRITT 3: Animationen erstellen (3 GIFs)"),
    ]
    
    success_count = 0
    
    for script, desc in scripts:
        if run_script(script, desc):
            success_count += 1
    
    # Zusammenfassung
    print("\n" + "="*70)
    print("  PIPELINE ABGESCHLOSSEN")
    print("="*70)
    print(f"\n  Erfolgreich: {success_count}/{len(scripts)} Schritte")
    
    # Zeige Ausgabedateien
    if success_count == len(scripts):
        print("\n  [OK] Alle Animationen erstellt:")
        
        files_to_check = [
            "output/01_orbit_classical.gif",
            "output/02_orbit_ssz_comparison.gif",
            "output/03_orbit_3d.gif",
            "orbit_data.json",
            "ssz_orbit_data.json"
        ]
        
        for file_path in files_to_check:
            full_path = Path(__file__).parent / file_path
            if full_path.exists():
                size = full_path.stat().st_size / 1024  # KB
                print(f"     • {file_path} ({size:.1f} KB)")
        
        print("\n  Öffne die GIFs im output/-Verzeichnis!")
    else:
        print("\n  [WARN] Einige Schritte sind fehlgeschlagen.")
        print("     Prüfe die Fehlermeldungen oben.")
    
    print("\n" + "="*70)


if __name__ == "__main__":
    main()
