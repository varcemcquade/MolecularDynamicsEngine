import sys, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

data = sys.argv[1] if len(sys.argv) > 1 else "energies.csv"
plots = Path(sys.argv[2] if len(sys.argv) > 2 else "plots")
plots.mkdir(exist_ok=True)

d = np.genfromtxt(data, delimiter=',', names=True)
TARGET_TEMP = 94.4

fig, axes = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
for ax, col, color in zip(axes, ["KE", "PE", "E_total"], ["#E69F00", "#56B4E9", "#009E73"]):
    ax.plot(d["step"], d[col], color=color, lw=0.6, alpha=0.85)
    ax.set_ylabel(col)
axes[-1].set_xlabel("Timestep")
axes[0].set_title("Energy vs. Time")
plt.tight_layout(); plt.savefig(plots / "energy_timeseries.png", dpi=150); plt.close()

m, b = np.polyfit(d["step"], d["E_total"], 1)
fit = m * d["step"] + b
r2 = 1 - np.sum((d["E_total"] - fit)**2) / np.sum((d["E_total"] - d["E_total"].mean())**2)
print(f"Drift: slope={m:+.4e} J/step  R²={r2:.6f}")

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(d["step"], d["E_total"], color="#009E73", lw=0.5, alpha=0.7)
ax.plot(d["step"], fit, color="#CC3311", lw=0.8, ls="--")
ax.set(title=f"Energy Drift  slope={m:+.3e}  R²={r2:.5f}", xlabel="Timestep", ylabel="E_total (J)")
plt.tight_layout(); plt.savefig(plots / "energy_drift.png", dpi=150); plt.close()

for col in ["KE", "PE", "E_total", "Temperature"]:
    print(f"  <{col}> = {d[col].mean():.6e} ± {d[col].std():.4e}")
T = d["Temperature"].mean()
print(f"  ΔT/T_target = {abs(T - TARGET_TEMP)/TARGET_TEMP*100:.2f}%")
