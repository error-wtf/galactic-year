"""
02_animate_orbit.py — Animierte GIF der Sonnen-Umlaufbahn in der Galaxie
=========================================================================

Liest die von 01_fetch_data.py erzeugten Parameter und:
  1. Integriert die 3D-Bahn der Sonne über mehrere galaktische Jahre
     mit scipy (Leapfrog + RK45) im MWPotential2014-ähnlichem Potential
  2. Rendert ein animiertes GIF mit matplotlib + Pillow

Galaktisches Potential (vereinfacht nach Bovy 2015 / MWPotential2014):
  Φ_total = Φ_bulge (Hernquist) + Φ_disk (Miyamoto-Nagai) + Φ_halo (NFW)

Ausgabe: output/solar_orbit.gif
         output/solar_orbit_3d.gif
"""

import json
import sys
import warnings
from pathlib import Path

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.animation import FuncAnimation, PillowWriter

warnings.filterwarnings("ignore")

# ── Pillow-Check ─────────────────────────────────────────────────────────────
try:
    from PIL import Image
    _has_pillow = True
except ImportError:
    _has_pillow = False
    print("WARNUNG: Pillow nicht gefunden → pip install Pillow")

# ── Scipy-Check ──────────────────────────────────────────────────────────────
try:
    from scipy.integrate import solve_ivp
    _has_scipy = True
except ImportError:
    _has_scipy = False
    print("WARNUNG: scipy nicht gefunden → Nutze eigenen Leapfrog-Integrator.")

# ── Pfade ────────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).parent
DATA_FILE  = BASE_DIR / "data" / "galactic_orbit_params.json"
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

GIF_2D  = OUTPUT_DIR / "solar_orbit_2d.gif"
GIF_3D  = OUTPUT_DIR / "solar_orbit_3d.gif"
PNG_NOW = OUTPUT_DIR / "solar_orbit_now.png"

# ── Styling ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor":  "#05060f",
    "axes.facecolor":    "#05060f",
    "axes.edgecolor":    "#2a2d4a",
    "text.color":        "#e0e4f0",
    "axes.labelcolor":   "#e0e4f0",
    "xtick.color":       "#8090b0",
    "ytick.color":       "#8090b0",
    "grid.color":        "#1a1d30",
    "grid.linewidth":    0.5,
    "font.family":       "monospace",
})

GOLD   = "#f5c842"
CYAN   = "#3ecfcf"
VIOLET = "#9b5de5"
WHITE  = "#e0e4f0"


# ═══════════════════════════════════════════════════════════════════════════
# Galaktisches Potential (vereinfacht nach MWPotential2014)
# Alle Größen in SI: m, kg, s
# ═══════════════════════════════════════════════════════════════════════════

KPC  = 3.085677581e19   # m
KM_S = 1.0e3            # m/s
G_SI = 6.67430e-11      # m³ kg⁻¹ s⁻²
M_SOL = 1.989e30        # kg
YR    = 365.25 * 24 * 3600   # s


def _hernquist_acc(pos, M, a):
    """Hernquist-Beschleunigung (Bulge). pos in m, M in kg, a in m."""
    r = np.linalg.norm(pos)
    return -G_SI * M / (r * (r + a)**2) * pos


def _miyamoto_nagai_acc(pos, M, a, b):
    """Miyamoto-Nagai-Beschleunigung (galaktische Disk). pos in m."""
    x, y, z = pos
    R2 = x**2 + y**2
    D  = np.sqrt(z**2 + b**2)
    denom = (R2 + (a + D)**2) ** 1.5
    axy   = -G_SI * M / denom
    az    = axy * (a + D) / D
    return np.array([axy * x, axy * y, az * z])


def _nfw_acc(pos, M200, r200, c):
    """NFW-Halo-Beschleunigung. c = Konzentrationsparameter."""
    rs  = r200 / c
    r   = np.linalg.norm(pos)
    rho0 = M200 / (4.0 * np.pi * rs**3 * (np.log(1 + c) - c / (1 + c)))
    M_r  = 4.0 * np.pi * rho0 * rs**3 * (np.log(1 + r / rs) - (r / rs) / (1 + r / rs))
    return -G_SI * M_r / r**2 * (pos / r)


def galactic_acceleration(pos, params):
    """Gesamte galaktische Beschleunigung an Position pos [m]."""
    M_bulge  = params["M_bulge"]
    a_bulge  = params["a_bulge"]
    M_disk   = params["M_disk"]
    a_disk   = params["a_disk"]
    b_disk   = params["b_disk"]
    M_halo   = params["M_halo"]
    r200     = params["r200"]
    c_halo   = params["c_halo"]

    acc  = _hernquist_acc(pos, M_bulge, a_bulge)
    acc += _miyamoto_nagai_acc(pos, M_disk, a_disk, b_disk)
    acc += _nfw_acc(pos, M_halo, r200, c_halo)
    return acc


# MWPotential2014-ähnliche Parameter (Bovy 2015)
MW_PARAMS = {
    "M_bulge":   3.0e10  * M_SOL,
    "a_bulge":   0.6     * KPC,
    "M_disk":    6.8e10  * M_SOL,
    "a_disk":    3.0     * KPC,
    "b_disk":    0.28    * KPC,
    "M_halo":    8.0e11  * M_SOL,
    "r200":      200.0   * KPC,
    "c_halo":    15.3,
}


# ═══════════════════════════════════════════════════════════════════════════
# Orbit-Integration
# ═══════════════════════════════════════════════════════════════════════════

def build_initial_conditions(full_data: dict) -> tuple:
    """
    Anfangsbedingungen der Sonne in galaktozentrischen Koordinaten (SI).
    Sonnenposition: (R0, 0, z_sun) [m]
    Sonnengeschwindigkeit: (U_sun, Theta0+V_sun, W_sun) [m/s]
    """
    # Neue JSON-Struktur: orbit_parameters (von 01_fetch_data.py)
    dp = full_data.get("orbit_parameters") or full_data.get("derived_parameters", {})
    R0    = dp["R0_kpc"] * KPC
    z0_pc = dp["vertical_oscillation"]["current_z_pc"]
    z0    = z0_pc * 3.085677581e16  # pc -> m

    Theta0 = dp["Theta0_km_s"] * KM_S
    U_sun  = dp["solar_peculiar_motion"]["U_km_s"] * KM_S
    V_sun  = dp["solar_peculiar_motion"]["V_km_s"] * KM_S
    W_sun  = dp["solar_peculiar_motion"]["W_km_s"] * KM_S

    pos = np.array([R0, 0.0, z0])
    vel = np.array([-U_sun, Theta0 + V_sun, W_sun])   # x→GC, y→Rotation, z→Nord
    return pos, vel


def leapfrog_integrate(pos0, vel0, dt_yr, n_steps):
    """
    Leapfrog-Integrator (symplektisch, energie­erhaltend).
    dt_yr   – Zeitschritt in Jahren
    n_steps – Anzahl Schritte
    """
    dt = dt_yr * YR
    pos  = pos0.copy()
    vel  = vel0.copy()
    traj = np.zeros((n_steps + 1, 6))
    traj[0, :3] = pos
    traj[0, 3:] = vel

    # Kick-Drift-Kick (Leapfrog)
    acc = galactic_acceleration(pos, MW_PARAMS)
    for i in range(n_steps):
        vel  = vel  + 0.5 * acc * dt
        pos  = pos  + vel * dt
        acc  = galactic_acceleration(pos, MW_PARAMS)
        vel  = vel  + 0.5 * acc * dt
        traj[i + 1, :3] = pos
        traj[i + 1, 3:] = vel

    return traj


def rk45_integrate(pos0, vel0, T_total_yr, n_points):
    """RK45-Integration via scipy (präziser, aber langsamer)."""
    def ode(t, y):
        pos = y[:3]
        vel = y[3:]
        acc = galactic_acceleration(pos, MW_PARAMS)
        return np.concatenate([vel, acc])

    t_span = (0, T_total_yr * YR)
    t_eval = np.linspace(0, T_total_yr * YR, n_points)
    y0     = np.concatenate([pos0, vel0])

    print(f"  scipy RK45: {T_total_yr / 1e6:.0f} Mio. Jahre, {n_points} Punkte ...")
    sol = solve_ivp(ode, t_span, y0, method="RK45", t_eval=t_eval,
                    rtol=1e-8, atol=1e-10, dense_output=False)
    traj = np.zeros((len(sol.t), 6))
    traj[:, :3] = sol.y[:3].T
    traj[:, 3:] = sol.y[3:].T
    return traj


# ═══════════════════════════════════════════════════════════════════════════
# Hilfsfunktionen: Spiralarm-Hintergrund
# ═══════════════════════════════════════════════════════════════════════════

def _draw_milky_way_background(ax):
    """Zeichnet schematische Spiralarme und galaktisches Zentrum."""
    # Hintergrund-Glow (galakt. Disk)
    for r, alpha in [(18, 0.04), (14, 0.06), (10, 0.09), (6, 0.12)]:
        circle = Circle((0, 0), r, color="#a0c0ff", alpha=alpha, transform=ax.transData)
        ax.add_patch(circle)

    # 4 Spiralarme (logarithmisch, vereinfacht)
    theta = np.linspace(0, 4 * np.pi, 500)
    for arm_offset, color, label in [
        (0.0,          "#4080ff", "Sgr Arm"),
        (np.pi / 2,    "#40c0ff", "Sct-Cen Arm"),
        (np.pi,        "#8040ff", "Perseus Arm"),
        (3 * np.pi / 2,"#40ffb0", "Lokaler Arm"),
    ]:
        k   = 0.22        # Spiralsteigung
        r0  = 2.5         # kpc Startradius
        r   = r0 * np.exp(k * theta)
        mask = r < 20
        theta_arm = theta[mask] + arm_offset
        r_arm     = r[mask]
        x_arm = r_arm * np.cos(theta_arm)
        y_arm = r_arm * np.sin(theta_arm)
        ax.plot(x_arm, y_arm, color=color, alpha=0.18, lw=1.8, zorder=1)

    # Galaktisches Zentrum
    ax.scatter([0], [0], s=350, color=GOLD, zorder=10,
               marker="*", edgecolors="#fff5a0", linewidths=0.5)
    ax.text(0, -1.3, "GC", color=GOLD, ha="center", va="top", fontsize=7, zorder=11)


# ═══════════════════════════════════════════════════════════════════════════
# 2D-Animation (galaktische Ebene)
# ═══════════════════════════════════════════════════════════════════════════

def create_2d_animation(traj_kpc, t_myr, params_meta):
    """Erstellt 2D-GIF der Sonnenbahn in der galaktischen Ebene (x-y)."""
    print("\n[GIF 1/2] Erstelle 2D-Animation ...")

    x = traj_kpc[:, 0]
    y = traj_kpc[:, 1]
    R0 = params_meta["R0_kpc"]
    n_frames = min(300, len(x))
    step     = max(1, len(x) // n_frames)
    indices  = np.arange(0, len(x), step)

    # Farbverlauf für die Spur (Zeit → Farbe)
    cmap = LinearSegmentedColormap.from_list(
        "orbit", ["#1a0050", VIOLET, CYAN, GOLD], N=256
    )

    fig, ax = plt.subplots(figsize=(8, 8), dpi=100,
                            facecolor="#05060f", subplot_kw={"facecolor": "#05060f"})
    ax.set_aspect("equal")
    lim = max(R0 * 1.6, 14.0)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(alpha=0.15, lw=0.4)
    ax.set_xlabel("x  [kpc]  (→ galakt. Zentrum)", color=WHITE, fontsize=8)
    ax.set_ylabel("y  [kpc]  (Rotationsrichtung)", color=WHITE, fontsize=8)
    ax.set_title("Sonnenbahn in der Milchstraße\n(galaktische Ebene)",
                 color=WHITE, fontsize=10, pad=10)

    _draw_milky_way_background(ax)

    # Vollständige Spur (blass)
    ax.plot(x, y, color=VIOLET, alpha=0.15, lw=0.7, zorder=2)

    # Aktive Spur (Farbverlauf)
    trail_len = max(20, n_frames // 8)
    trail_line, = ax.plot([], [], lw=1.5, color=CYAN, alpha=0.7, zorder=3)
    sun_dot     = ax.scatter([], [], s=120, color=GOLD, zorder=12,
                              edgecolors="white", linewidths=0.8)
    time_text   = ax.text(0.02, 0.96, "", transform=ax.transAxes,
                           color=WHITE, fontsize=8, va="top",
                           bbox=dict(boxstyle="round,pad=0.3", fc="#0a0c1a", alpha=0.7))

    # Aktuelle Sonnenposition markieren
    ax.scatter([x[0]], [y[0]], s=80, color=GOLD, marker="o",  # ⊙
               zorder=13, edgecolors="white", linewidths=0.5)

    R_circ = Circle((0, 0), R0, color=CYAN, fill=False,
                    lw=0.6, alpha=0.3, linestyle="--", zorder=2)
    ax.add_patch(R_circ)
    ax.text(R0 * 0.72, R0 * 0.72, f"R₀={R0:.2f} kpc",
            color=CYAN, fontsize=6.5, alpha=0.6, zorder=3)

    # Parameterbox
    T_gyr  = t_myr[-1] / 1e3
    box_txt = (
        f"R₀   = {R0:.3f} kpc  (GRAVITY 2019)\n"
        f"Θ₀   = {params_meta['Theta0_km_s']:.0f} km/s  (Reid+2014)\n"
        f"T_gal= {params_meta['T_gal_myr']:.0f} Mio. a\n"
        f"Dauer= {T_gyr:.2f} Gyr"
    )
    ax.text(0.97, 0.03, box_txt, transform=ax.transAxes, color=WHITE,
            fontsize=6.5, ha="right", va="bottom", family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", fc="#0a0c1a", alpha=0.8))

    def init():
        trail_line.set_data([], [])
        sun_dot.set_offsets(np.empty((0, 2)))
        time_text.set_text("")
        return trail_line, sun_dot, time_text

    def update(frame):
        idx  = indices[frame]
        i0   = max(0, idx - trail_len * step)
        i1   = idx + 1
        xtr  = x[i0:i1:step]
        ytr  = y[i0:i1:step]
        trail_line.set_data(xtr, ytr)
        sun_dot.set_offsets([[x[idx], y[idx]]])
        age_myr = t_myr[idx]
        time_text.set_text(
            f"t = {age_myr / 1e3:.3f} Gyr\n"
            f"    ({age_myr:.0f} Mio. Jahre)"
        )
        return trail_line, sun_dot, time_text

    anim = FuncAnimation(fig, update, frames=len(indices),
                          init_func=init, blit=True, interval=30)

    writer = PillowWriter(fps=30)
    anim.save(str(GIF_2D), writer=writer,
              savefig_kwargs={"facecolor": "#05060f"})
    plt.close(fig)
    print(f"  [OK] {GIF_2D}")


# ═══════════════════════════════════════════════════════════════════════════
# 3D-Animation (x-y-z, zeigt vertikale Oszillation)
# ═══════════════════════════════════════════════════════════════════════════

def create_3d_animation(traj_kpc, t_myr, params_meta):
    """Erstellt 3D-GIF mit galktischer Ebene + vertikaler Oszillation."""
    print("\n[GIF 2/2] Erstelle 3D-Animation ...")
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    x = traj_kpc[:, 0]
    y = traj_kpc[:, 1]
    z = traj_kpc[:, 2]

    n_frames = min(200, len(x))
    step     = max(1, len(x) // n_frames)
    indices  = np.arange(0, len(x), step)

    R0  = params_meta["R0_kpc"]
    lim = max(R0 * 1.5, 12.0)
    zlim = max(abs(z).max() * 1.5, 0.3)

    fig = plt.figure(figsize=(9, 8), facecolor="#05060f")
    ax  = fig.add_subplot(111, projection="3d", facecolor="#05060f")
    ax.set_facecolor("#05060f")

    # Galaktische Ebene (transparente Disk)
    theta_tmp = np.linspace(0, 2 * np.pi, 60)
    r_tmp     = np.linspace(0, lim, 10)
    T2, R2    = np.meshgrid(theta_tmp, r_tmp)
    Xd = R2 * np.cos(T2)
    Yd = R2 * np.sin(T2)
    Zd = np.zeros_like(Xd)
    ax.plot_surface(Xd, Yd, Zd, alpha=0.06, color="#4060c0", zorder=0)

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-zlim, zlim)
    ax.set_xlabel("x [kpc]", color=WHITE, fontsize=7)
    ax.set_ylabel("y [kpc]", color=WHITE, fontsize=7)
    ax.set_zlabel("z [kpc]", color=WHITE, fontsize=7)
    ax.set_title("Sonnen-Orbit 3D\n(inkl. vertikale Oszillation)",
                 color=WHITE, fontsize=9)
    ax.tick_params(colors=WHITE, labelsize=6)

    # Vollständige Spur (blass)
    ax.plot(x, y, z, color=VIOLET, alpha=0.12, lw=0.6)

    trail_len = max(30, n_frames // 6)
    trail_line, = ax.plot([], [], [], lw=1.8, color=CYAN, alpha=0.8)
    sun_dot     = ax.scatter([], [], [], s=150, color=GOLD, zorder=15,
                              edgecolors="white", linewidths=0.8)
    time_text   = ax.text2D(0.02, 0.95, "", transform=ax.transAxes,
                             color=WHITE, fontsize=7)

    # GC-Stern
    ax.scatter([0], [0], [0], s=300, color=GOLD, marker="*", zorder=20)

    def init():
        trail_line.set_data_3d([], [], [])
        sun_dot._offsets3d = (np.array([]), np.array([]), np.array([]))
        time_text.set_text("")
        return trail_line, sun_dot, time_text

    def update(frame):
        idx = indices[frame]
        i0  = max(0, idx - trail_len * step)
        i1  = idx + 1
        xtr = x[i0:i1:step]
        ytr = y[i0:i1:step]
        ztr = z[i0:i1:step]
        trail_line.set_data_3d(xtr, ytr, ztr)
        sun_dot._offsets3d = (np.array([x[idx]]),
                               np.array([y[idx]]),
                               np.array([z[idx]]))
        time_text.set_text(f"t = {t_myr[idx] / 1e3:.3f} Gyr  "
                            f"z = {z[idx] * 1000:.1f} pc")
        # Kamera leicht rotieren
        ax.view_init(elev=25, azim=frame * 360 / n_frames)
        return trail_line, sun_dot, time_text

    anim = FuncAnimation(fig, update, frames=len(indices),
                          init_func=init, blit=False, interval=40)
    writer = PillowWriter(fps=25)
    anim.save(str(GIF_3D), writer=writer,
              savefig_kwargs={"facecolor": "#05060f"})
    plt.close(fig)
    print(f"  [OK] {GIF_3D}")


# ═══════════════════════════════════════════════════════════════════════════
# Statisches Überblicks-PNG
# ═══════════════════════════════════════════════════════════════════════════

def create_overview_png(traj_kpc, t_myr, params_meta):
    """4-Panel Übersichtsgrafik: Bahn, r(t), z(t), v(t)."""
    print("\n[PNG]  Erstelle Uebersichts-PNG ...")
    x = traj_kpc[:, 0]
    y = traj_kpc[:, 1]
    z = traj_kpc[:, 2]
    vx = traj_kpc[:, 3] if traj_kpc.shape[1] > 3 else np.zeros_like(x)

    r = np.sqrt(x**2 + y**2)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10),
                              facecolor="#05060f",
                              gridspec_kw={"hspace": 0.35, "wspace": 0.35})

    # Panel 1: Bahn in x-y
    ax = axes[0, 0]
    ax.set_facecolor("#05060f")
    _draw_milky_way_background(ax)
    sc = ax.scatter(x[::5], y[::5],
                     c=t_myr[::5], cmap="plasma", s=1.5, alpha=0.6, zorder=5)
    R0 = params_meta["R0_kpc"]
    ax.add_patch(Circle((0, 0), R0, color=CYAN, fill=False, lw=0.7, ls="--", alpha=0.4))
    ax.scatter([x[0]], [y[0]], s=120, color=GOLD, zorder=15, marker="*",
               edgecolors="white", linewidths=0.7, label="Heute")
    ax.set_aspect("equal")
    lim = max(R0 * 1.7, 14)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("x [kpc]", color=WHITE); ax.set_ylabel("y [kpc]", color=WHITE)
    ax.set_title("Bahn in der galaktischen Ebene", color=WHITE, fontsize=9)
    plt.colorbar(sc, ax=ax, label="t [Mio. Jahre]", fraction=0.04)
    ax.legend(fontsize=7, facecolor="#1a1d30", labelcolor=WHITE)

    # Panel 2: r(t)
    ax = axes[0, 1]
    ax.set_facecolor("#05060f")
    ax.plot(t_myr / 1e3, r, color=CYAN, lw=0.9, alpha=0.9)
    ax.axhline(R0, color=GOLD, ls="--", lw=0.8, alpha=0.6, label=f"R₀={R0:.2f} kpc")
    ax.set_xlabel("t [Gyr]", color=WHITE); ax.set_ylabel("R  [kpc]", color=WHITE)
    ax.set_title("Galaktozentrischer Abstand R(t)", color=WHITE, fontsize=9)
    ax.legend(fontsize=7, facecolor="#1a1d30", labelcolor=WHITE)
    ax.grid(alpha=0.2)

    # Panel 3: z(t) – vertikale Oszillation
    ax = axes[1, 0]
    ax.set_facecolor("#05060f")
    ax.plot(t_myr / 1e3, z * 1000, color=VIOLET, lw=0.9, alpha=0.9)
    ax.axhline(0, color=WHITE, ls="--", lw=0.5, alpha=0.4, label="Galakt. Ebene")
    ax.set_xlabel("t [Gyr]", color=WHITE)
    ax.set_ylabel("z  [pc]", color=WHITE)
    ax.set_title("Vertikale Oszillation z(t)", color=WHITE, fontsize=9)
    ax.legend(fontsize=7, facecolor="#1a1d30", labelcolor=WHITE)
    ax.grid(alpha=0.2)

    # Panel 4: Rotationskurve (aus Daten + Sonnenposition)
    ax = axes[1, 1]
    ax.set_facecolor("#05060f")
    # Analytische Kurve
    r_arr = np.linspace(1, 25, 300)
    v_arr = 240.0 - 0.2 * (r_arr - 8.34)  # leicht fallend, nach Reid+2014
    ax.plot(r_arr, v_arr, color=CYAN, lw=1.5, label="Rotationskurve (Reid+2014)")
    ax.scatter([R0], [params_meta["Theta0_km_s"]], s=120, color=GOLD,
               zorder=10, marker="o", edgecolors="white", lw=0.7, label="Sonne (☉)")
    ax.set_xlabel("R  [kpc]", color=WHITE)
    ax.set_ylabel("v_circ  [km/s]", color=WHITE)
    ax.set_title("Galaktische Rotationskurve", color=WHITE, fontsize=9)
    ax.legend(fontsize=7, facecolor="#1a1d30", labelcolor=WHITE)
    ax.set_xlim(0, 26); ax.set_ylim(180, 280)
    ax.grid(alpha=0.2)

    # Haupttitel
    T_gyr = t_myr[-1] / 1e3
    n_orb  = T_gyr * 1e9 / (params_meta["T_gal_myr"] * 1e6)
    fig.suptitle(
        f"Sonnenbahn in der Milchstraße  |  {T_gyr:.1f} Gyr  ≈  {n_orb:.1f} Umläufe\n"
        f"R₀={R0:.3f} kpc  Θ₀={params_meta['Theta0_km_s']:.0f} km/s  "
        f"T_gal={params_meta['T_gal_myr']:.0f} Mio. Jahre",
        color=WHITE, fontsize=11, y=0.995,
    )

    fig.savefig(str(PNG_NOW), dpi=150, facecolor="#05060f", bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {PNG_NOW}")


# ═══════════════════════════════════════════════════════════════════════════
# Hauptprogramm
# ═══════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 65)
    print("  GALAXY-YEAR | Orbit-Animation")
    print("=" * 65)

    # ── Parameter laden ──────────────────────────────────────────────────
    if DATA_FILE.exists():
        print(f"\nLade Parameter aus {DATA_FILE} …")
        with open(DATA_FILE, encoding="utf-8") as f:
            data = json.load(f)
        dp  = data["orbit_parameters"]
        # Kompatibilitaet: Anfangsbedingungen aus data direkt
        data["derived_parameters"] = dp  # Alias fuer build_initial_conditions()
    else:
        print(f"\n[!] {DATA_FILE} nicht gefunden.")
        print("   Bitte zuerst `python 01_fetch_data.py` ausfuehren!")
        print("   Nutze eingebettete Standardwerte (Reid+2014 / GRAVITY 2019) ...\n")
        dp = {
            "R0_kpc": 8.122, "Theta0_km_s": 240.0,
            "galactic_year_myr": 212.0,
            "solar_peculiar_motion": {"U_km_s": 11.1, "V_km_s": 12.24, "W_km_s": 7.25},
            "vertical_oscillation": {"amplitude_kpc": 0.093, "period_myr": 70, "current_z_pc": 17},
        }
        data = {"orbit_parameters": dp, "derived_parameters": dp}

    params_meta = {
        "R0_kpc":       dp["R0_kpc"],
        "Theta0_km_s":  dp["Theta0_km_s"],
        "T_gal_myr":    dp.get("galactic_year_myr", 212.0),
    }

    # ── Anfangsbedingungen ───────────────────────────────────────────────
    print("\nBerechne Anfangsbedingungen ...")
    pos0, vel0 = build_initial_conditions(data)
    print(f"  pos0 = ({pos0[0]/KPC:.3f}, {pos0[1]/KPC:.3f}, {pos0[2]/KPC*1000:.1f} pc) kpc")
    print(f"  vel0 = ({vel0[0]/KM_S:.1f}, {vel0[1]/KM_S:.1f}, {vel0[2]/KM_S:.1f}) km/s")

    # ── Orbit integrieren ────────────────────────────────────────────────
    # 5 Gyr = ~23 galaktische Umläufe
    T_total_yr = 5.0e9
    n_points   = 8000

    if _has_scipy:
        print(f"\nIntegriere Orbit ({T_total_yr/1e9:.0f} Gyr, {n_points} Punkte) ...")
        traj = rk45_integrate(pos0, vel0, T_total_yr, n_points)
    else:
        print(f"\nLeapfrog-Integration ({T_total_yr/1e9:.0f} Gyr, {n_points} Schritte) ...")
        dt_yr = T_total_yr / n_points
        traj  = leapfrog_integrate(pos0, vel0, dt_yr, n_points)

    # ── Einheiten → kpc ──────────────────────────────────────────────────
    traj_kpc = traj.copy()
    traj_kpc[:, :3] /= KPC      # m → kpc
    traj_kpc[:, 3:] /= KM_S     # m/s → km/s

    t_myr = np.linspace(0, T_total_yr / 1e6 / YR * YR, n_points + 1) / (YR * 1e6)
    # Korrekte Zeitachse
    t_myr = np.linspace(0, T_total_yr / YR / 1e6, len(traj_kpc))

    R = np.sqrt(traj_kpc[:, 0]**2 + traj_kpc[:, 1]**2)
    print("\n  Orbit-Statistik:")
    print(f"  R_min = {R.min():.3f} kpc   R_max = {R.max():.3f} kpc")
    print(f"  z_max = {np.abs(traj_kpc[:, 2]).max() * 1000:.1f} pc")

    # ── Grafiken erstellen ────────────────────────────────────────────────
    create_overview_png(traj_kpc, t_myr, params_meta)
    create_2d_animation(traj_kpc, t_myr, params_meta)
    create_3d_animation(traj_kpc, t_myr, params_meta)

    print(f"\n{'=' * 65}")
    print(f"  [OK] Fertig!")
    print(f"  -> {GIF_2D}")
    print(f"  -> {GIF_3D}")
    print(f"  -> {PNG_NOW}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
