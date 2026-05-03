import sys, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

data = sys.argv[1] if len(sys.argv) > 1 else "msd.csv"
plots = Path(sys.argv[2] if len(sys.argv) > 2 else "plots")
plots.mkdir(exist_ok=True)

raw = np.genfromtxt(data, delimiter=',', skip_header=1)
t = raw[:, 0] * 1e-15      # s
msd = raw[:, 1] * 1e-20    # m²
msd_A2 = raw[:, 1]         # Å²

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t * 1e12, msd_A2, color="#0072B2", lw=0.8)
ax.set(title="MSD vs. Time", xlabel="Time (ps)", ylabel="MSD (Å²)")
plt.tight_layout(); plt.savefig(plots / "msd_timeseries.png", dpi=150); plt.close()

lin_t, lin_msd = t[len(t)//5:], msd[len(msd)//5:]
m, b = np.polyfit(lin_t, lin_msd, 1)
D = m / 6
fit = m * t + b
r2 = 1 - np.sum((lin_msd - (m*lin_t+b))**2) / np.sum((lin_msd - lin_msd.mean())**2)
print(f"D = {D:.4e} m²/s  R²={r2:.6f}  vs exp ratio={D/2.05e-9:.3f}")

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(t * 1e12, msd * 1e20, color="#0072B2", lw=0.8, alpha=0.7)
ax.plot(t * 1e12, fit * 1e20, color="#CC3311", lw=0.9, ls="--")
ax.annotate(f"D={D:.3e} m²/s\nR²={r2:.4f}", xy=(0.05, 0.85), xycoords="axes fraction", color="#CC3311")
ax.set(title="MSD — Einstein Fit", xlabel="Time (ps)", ylabel="MSD (Å²)")
plt.tight_layout(); plt.savefig(plots / "msd_einstein_fit.png", dpi=150); plt.close()

pos = (msd > 0) & (t > 0)
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(np.log10(t[pos]), np.log10(msd[pos]), color="#0072B2", lw=0.8)
ax.set(title="MSD — Log-Log", xlabel="log₁₀(t / s)", ylabel="log₁₀(MSD / m²)")
plt.tight_layout(); plt.savefig(plots / "msd_loglog.png", dpi=150); plt.close()
