import pathlib
txt = pathlib.Path("02_animate_orbit.py").read_text(encoding="utf-8")
bad = "\u2609"   # ☉ sun symbol
txt = txt.replace(f'marker="{bad}"', 'marker="o"')
pathlib.Path("02_animate_orbit.py").write_text(txt, encoding="utf-8")
print("Fixed", txt.count('marker="o"'), "markers")
