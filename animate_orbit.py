"""
animate_orbit.py

Erstellt animierte GIFs der Umlaufbahn:
1. Klassische Kepler-Bahn
2. SSZ-korrigierte Bahn (Segmented Spacetime)

Beide Animationen zeigen den Vergleich und die Unterschiede.
"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import imageio
from pathlib import Path
import json
from tqdm import tqdm
from typing import Dict, Optional


def load_data(filename: str) -> Optional[Dict]:
    """Lädt Orbitdaten aus JSON."""
    data_path = Path(__file__).parent / filename
    if not data_path.exists():
        return None
    with open(data_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def create_classical_animation(data: Dict, output_filename: str = "orbit_classical.gif",
                                  fps: int = 30, duration: int = 6) -> Path:
    """
    Erstellt Animation der klassischen Kepler-Bahn.
    """
    x = np.array(data['data']['x'])
    y = np.array(data['data']['y'])
    z = np.array(data['data']['z'])
    t = np.array(data['data']['t'])
    params = data['metadata']['parameters']
    
    n_frames = fps * duration
    frame_indices = np.linspace(0, len(x) - 1, n_frames, dtype=int)
    
    fig, ax = plt.subplots(figsize=(11, 10))
    
    # Elemente
    full_orbit, = ax.plot([], [], 'b-', alpha=0.25, linewidth=1.5, label='Umlaufbahn')
    completed, = ax.plot([], [], 'b-', linewidth=2.5, alpha=0.9)
    sun, = ax.plot([], [], 'yo', markersize=14, markeredgecolor='orange',
                   markeredgewidth=2.5, zorder=5)
    center, = ax.plot([0], [0], 'k*', markersize=22, label='Sgr A*', zorder=5)
    
    # Peri- und Apogalaktikum markieren
    r = np.sqrt(x**2 + y**2)
    peri_idx = np.argmin(r)
    apo_idx = np.argmax(r)
    ax.plot(x[peri_idx], y[peri_idx], 'rv', markersize=10, label='Perigalaktikum')
    ax.plot(x[apo_idx], y[apo_idx], 'r^', markersize=10, label='Apogalaktikum')
    
    # Text
    time_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, fontsize=11,
                        verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    info_text = ax.text(0.98, 0.02, '', transform=ax.transAxes, fontsize=10,
                        verticalalignment='bottom', horizontalalignment='right',
                        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
    
    # Layout
    max_range = max(np.max(np.abs(x)), np.max(np.abs(y))) * 1.15
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_aspect('equal')
    ax.set_xlabel('x [kpc]', fontsize=12)
    ax.set_ylabel('y [kpc]', fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='lower right', fontsize=10, framealpha=0.9)
    
    ax.set_title('Klassische Kepler-Umlaufbahn der Sonne\n' +
                 f'R₀={params["R_0"]:.1f} kpc, e={params["e"]:.3f}, T={params["T_orbit"]/1e6:.0f} Myr',
                 fontsize=14, fontweight='bold')
    
    # Mittlere Bahn
    circle = plt.Circle((0, 0), params["R_0"], fill=False, linestyle='--',
                        color='gray', alpha=0.4, linewidth=1.5)
    ax.add_patch(circle)
    
    frames = []
    print(f"Erstelle klassische Animation ({n_frames} Frames)...")
    
    for frame_idx in tqdm(frame_indices):
        full_orbit.set_data(x, y)
        completed.set_data(x[:frame_idx+1], y[:frame_idx+1])
        sun.set_data([x[frame_idx]], [y[frame_idx]])
        
        t_my = t[frame_idx] / 1e6
        r_kpc = np.sqrt(x[frame_idx]**2 + y[frame_idx]**2)
        
        time_text.set_text(f'Zeit: {t_my:.1f} Myr')
        info_text.set_text(f'r = {r_kpc:.2f} kpc\n' +
                          f'z = {z[frame_idx]:.3f} kpc')
        
        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))[:,:,:3]
        frames.append(frame)
    
    plt.close(fig)
    
    output_path = Path(__file__).parent / "output" / output_filename
    output_path.parent.mkdir(exist_ok=True)
    imageio.mimsave(output_path, frames, fps=fps)
    print(f"[OK] Gespeichert: {output_path}")
    
    return output_path


def create_ssz_animation(classical_data: Dict, ssz_data: Dict,
                          output_filename: str = "orbit_ssz.gif",
                          fps: int = 30, duration: int = 6) -> Path:
    """
    Erstellt Animation mit SSZ-Vergleich.
    Zeigt klassische vs SSZ-Bahn nebeneinander.
    """
    # Klassische Daten
    x_cl = np.array(classical_data['data']['x'])
    y_cl = np.array(classical_data['data']['y'])
    z_cl = np.array(classical_data['data']['z'])
    t_cl = np.array(classical_data['data']['t'])
    
    # SSZ Daten
    x_ssz = np.array(ssz_data['data']['x'])
    y_ssz = np.array(ssz_data['data']['y'])
    z_ssz = np.array(ssz_data['data']['z'])
    xi = np.array(ssz_data['data']['xi'])
    D = np.array(ssz_data['data']['D'])
    
    params = classical_data['metadata']['parameters']
    comp = ssz_data['metadata']['comparison']
    
    n_frames = fps * duration
    frame_indices = np.linspace(0, len(x_cl) - 1, n_frames, dtype=int)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    max_range = max(np.max(np.abs(x_cl)), np.max(np.abs(y_cl))) * 1.15
    
    for ax in [ax1, ax2]:
        ax.set_xlim(-max_range, max_range)
        ax.set_ylim(-max_range, max_range)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.set_xlabel('x [kpc]', fontsize=11)
        ax.set_ylabel('y [kpc]', fontsize=11)
        center, = ax.plot([0], [0], 'k*', markersize=18, zorder=5)
        circle = plt.Circle((0, 0), params["R_0"], fill=False, linestyle='--',
                           color='gray', alpha=0.4)
        ax.add_patch(circle)
    
    ax1.set_title('Klassische Kepler-Bahn', fontsize=13, fontweight='bold')
    ax2.set_title('SSZ-Korrigierte Bahn', fontsize=13, fontweight='bold', color='darkgreen')
    
    # Klassische Elemente
    full_cl, = ax1.plot([], [], 'b-', alpha=0.25, linewidth=1.5)
    comp_cl, = ax1.plot([], [], 'b-', linewidth=2.5, alpha=0.9)
    sun_cl, = ax1.plot([], [], 'yo', markersize=12, markeredgecolor='orange', markeredgewidth=2)
    
    # SSZ Elemente
    full_ssz, = ax2.plot([], [], 'g-', alpha=0.25, linewidth=1.5)
    comp_ssz, = ax2.plot([], [], 'g-', linewidth=2.5, alpha=0.9)
    sun_ssz, = ax2.plot([], [], 'go', markersize=12, markeredgecolor='darkgreen', markeredgewidth=2)
    
    # SSZ Info
    xi_text = ax2.text(0.98, 0.98, '', transform=ax2.transAxes, fontsize=10,
                        verticalalignment='top', horizontalalignment='right',
                        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))
    
    # Zeit
    time_text = fig.text(0.5, 0.95, '', ha='center', fontsize=14, fontweight='bold',
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    fig.suptitle(f'Sonnenumlaufbahn: Klassisch vs SSZ\n' +
                 f'R₀={params["R_0"]:.1f} kpc, φ={ssz_data["metadata"]["ssz_parameters"]["phi"]:.6f}',
                 fontsize=15, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    frames = []
    print(f"Erstelle SSZ-Vergleichs-Animation ({n_frames} Frames)...")
    
    for frame_idx in tqdm(frame_indices):
        # Klassische Bahn
        full_cl.set_data(x_cl, y_cl)
        comp_cl.set_data(x_cl[:frame_idx+1], y_cl[:frame_idx+1])
        sun_cl.set_data([x_cl[frame_idx]], [y_cl[frame_idx]])
        
        # SSZ Bahn
        full_ssz.set_data(x_ssz, y_ssz)
        comp_ssz.set_data(x_ssz[:frame_idx+1], y_ssz[:frame_idx+1])
        sun_ssz.set_data([x_ssz[frame_idx]], [y_ssz[frame_idx]])
        
        t_my = t_cl[frame_idx] / 1e6
        time_text.set_text(f'Zeit: {t_my:.1f} Millionen Jahre')
        
        xi_text.set_text(f'Ξ = {xi[frame_idx]:.6f}\n' +
                        f'D = {D[frame_idx]:.6f}\n' +
                        f'z = {z_ssz[frame_idx]:.3f} kpc')
        
        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))[:,:,:3]
        frames.append(frame)
    
    plt.close(fig)
    
    output_path = Path(__file__).parent / "output" / output_filename
    output_path.parent.mkdir(exist_ok=True)
    imageio.mimsave(output_path, frames, fps=fps)
    print(f"[OK] Gespeichert: {output_path}")
    
    return output_path


def create_3d_animation(data: Dict, output_filename: str = "orbit_3d.gif",
                        fps: int = 30, duration: int = 6) -> Path:
    """
    Erstellt 3D-Animation mit rotierender Ansicht.
    """
    from mpl_toolkits.mplot3d import Axes3D
    
    x = np.array(data['data']['x'])
    y = np.array(data['data']['y'])
    z = np.array(data['data']['z'])
    t = np.array(data['data']['t'])
    params = data['metadata']['parameters']
    
    n_frames = fps * duration
    frame_indices = np.linspace(0, len(x) - 1, n_frames, dtype=int)
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Elemente
    full_line, = ax.plot([], [], [], 'b-', alpha=0.3, linewidth=1)
    comp_line, = ax.plot([], [], [], 'b-', linewidth=2)
    sun_dot, = ax.plot([], [], [], 'yo', markersize=10, markeredgecolor='orange')
    center_dot, = ax.plot([0], [0], [0], 'k*', markersize=15)
    
    # Galaktische Ebene
    max_range = max(np.max(np.abs(x)), np.max(np.abs(y))) * 1.1
    xx, yy = np.meshgrid(np.linspace(-max_range, max_range, 10),
                          np.linspace(-max_range, max_range, 10))
    zz = np.zeros_like(xx)
    ax.plot_surface(xx, yy, zz, alpha=0.1, color='gray')
    
    ax.set_xlabel('x [kpc]')
    ax.set_ylabel('y [kpc]')
    ax.set_zlabel('z [kpc]')
    ax.set_title('3D-Umlaufbahn mit vertikaler Oszillation', fontsize=13, fontweight='bold')
    
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-0.3, 0.3)
    
    frames = []
    print(f"Erstelle 3D-Animation ({n_frames} Frames)...")
    
    for i, frame_idx in enumerate(tqdm(frame_indices)):
        full_line.set_data(x, y)
        full_line.set_3d_properties(z)
        
        comp_line.set_data(x[:frame_idx+1], y[:frame_idx+1])
        comp_line.set_3d_properties(z[:frame_idx+1])
        
        sun_dot.set_data([x[frame_idx]], [y[frame_idx]])
        sun_dot.set_3d_properties([z[frame_idx]])
        
        # Rotierende Ansicht
        angle = 30 + (i / n_frames) * 360
        ax.view_init(elev=25, azim=angle)
        
        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))[:,:,:3]
        frames.append(frame)
    
    plt.close(fig)
    
    output_path = Path(__file__).parent / "output" / output_filename
    output_path.parent.mkdir(exist_ok=True)
    imageio.mimsave(output_path, frames, fps=fps)
    print(f"[OK] Gespeichert: {output_path}")
    
    return output_path


def main():
    print("="*70)
    print("ORBIT-ANIMATION ERSTELLEN")
    print("="*70)
    
    # Lade Daten
    classical_data = load_data("orbit_data.json")
    ssz_data = load_data("ssz_orbit_data.json")
    
    if classical_data is None:
        print("\n[WARN] Keine klassischen Orbitdaten gefunden!")
        print("  Führe zuerst aus: python fetch_orbit_data.py")
        return
    
    if ssz_data is None:
        print("\n[WARN] Keine SSZ-Daten gefunden!")
        print("  Führe zuerst aus: python calculate_ssz_orbit.py")
        return
    
    print(f"\n[OK] Daten geladen:")
    print(f"  Klassisch: {classical_data['metadata']['n_points']} Punkte")
    print(f"  SSZ:       {ssz_data['metadata']['n_points']} Punkte")
    
    # Erstelle Animationen
    print("\n" + "="*70)
    
    gif1 = create_classical_animation(classical_data, 
                                       output_filename="01_orbit_classical.gif")
    
    gif2 = create_ssz_animation(classical_data, ssz_data,
                                 output_filename="02_orbit_ssz_comparison.gif")
    
    gif3 = create_3d_animation(classical_data,
                              output_filename="03_orbit_3d.gif")
    
    print("\n" + "="*70)
    print("ALLE ANIMATIONEN ERSTELLT!")
    print("="*70)
    print(f"\nAusgabedateien:")
    print(f"  1. {gif1.name}")
    print(f"  2. {gif2.name}")
    print(f"  3. {gif3.name}")
    print(f"\nVerzeichnis: {gif1.parent}")


if __name__ == "__main__":
    main()
