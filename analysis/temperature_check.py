import sys, numpy as np, matplotlib.pyplot as plt
from pathlib import Path

data = sys.argv[1] if len(sys.argv) > 1 else "energies.csv"
plots = Path(sys.argv[2] if len(sys.argv) > 2 else "plots")
plots.mkdir(exist_ok=True)

d = np.genfromtxt(data, delimiter=',', names=True)
T, TARGET, WINDOW = d["Temperature"], 94.4, 20
roll = np.convolve(T, np.ones(WINDOW)/WINDOW, mode='same')

fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(d["step"], T,    color="#999", lw=0.4, alpha=0.6)
ax.plot(d["step"], roll, color="#E69F00", lw=0.9)
ax.axhline(TARGET, color="#CC3311", lw=0.8, ls="--", label=f"Target {TARGET} K")
ax.set(title="Temperature vs. Time", xlabel="Timestep", ylabel="Temperature (K)")
ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(plots / "temperature_timeseries.png", dpi=150); plt.close()

print(f"mean={T.mean():.4f}  SD={T.std():.4f}  ΔT={abs(T.mean()-TARGET):.4f} K")

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(T, bins=40, density=True, color="#56B4E9", alpha=0.8, edgecolor="white")
xd = np.linspace(T.min(), T.max(), 200)
ax.plot(xd, np.exp(-0.5*((xd-T.mean())/T.std())**2)/(T.std()*np.sqrt(2*np.pi)), color="#0072B2", lw=1)
ax.axvline(TARGET,   color="#CC3311", lw=0.9, ls="--", label="Target")
ax.axvline(T.mean(), color="#009E73", lw=0.9,           label="Mean")
ax.set(title=f"Temperature Distribution  mean={T.mean():.2f} K", xlabel="K", ylabel="Density")
ax.legend(fontsize=9)
plt.tight_layout(); plt.savefig(plots / "temperature_distribution.png", dpi=150); plt.close()

n = len(T)//2
print(f"Equilibration: first={T[:n].mean():.4f}  second={T[n:].mean():.4f} K")
