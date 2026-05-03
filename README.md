######## MolecularDynamicsEngine ########

Simple MD simulation engine using Velocity-Verlet integrator in C++/CUDA with Python analysis scripts.

Requirements:

- MSVC (Visual Studio 2022+)
- CUDA Toolkit 12.6
- CMake 3.18+
- Python 3 — `py -m pip install -r requirements.txt`

Build:

```bash
cmake --preset x64-debug
cmake --build out/build/x64-debug
```

Or use **Ctrl+Shift+B** in VS Code.

Run:

```bash
./out/build/x64-debug/MolecularDynamicsEngine.exe
```

Outputs: `energies.csv`, `msd.csv`, `trajectory.dcd`

Analysis:

```bash
# Run all scripts
py analysis/analysis.py

# With custom paths (positional: energies msd plots)
py analysis/analysis.py energies.csv msd.csv plots/
```

Plots saved to `analysis/plots/`.