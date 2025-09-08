import pandas as pd
from scipy.stats import qmc

# import numpy as np


param_ranges = {
    # "rho": (7.8e-9, 7.8e-9),  # Density [Tonne/mm^3]
    "E": (180e3, 220e3),  # Young's modulus [MPa]
    "nu": (0.25, 0.35),  # Poisson's ratio [-]
    "A": (400, 700),  # JC hardening A (yield strength) [MPa]
    "B": (800, 1200),  # JC hardening B [MPa]
    "n": (0.15, 0.35),  # JC hardening exponent [-]
    "D1": (1, 1.2),  # JC damage parameter [-]
    "D2": (0.05, 0.15),  # JC damage parameter [-]
    "D3": (-0.7, -0.3),  # JC damage parameter [-]
    "uts": (1000, 1300),  # Ultimate tensile strength [MPa]
    "kc": (400, 6000),  # Fracture toughness [MPa mm^1/2]
}

n_samples = 1
output_file_csv = "material_parameters.csv"
output_file_py = "material_parameters.py"


dim = len(param_ranges)
sampler = qmc.LatinHypercube(d=dim, seed=42)
sample = sampler.random(n=n_samples)

l_bounds = [low for low, high in param_ranges.values()]
u_bounds = [high for low, high in param_ranges.values()]
scaled_sample = qmc.scale(sample, l_bounds, u_bounds)


df = pd.DataFrame(scaled_sample, columns=param_ranges.keys())
df.insert(0, "id", [f"{i:04d}" for i in range(1, n_samples + 1)])
df.insert(1, "rho", 7.8e-9)


df = df.round(
    {
        "E": 0,
        "nu": 3,
        "A": 0,
        "B": 0,
        "n": 3,
        "D1": 3,
        "D2": 3,
        "D3": 3,
        "uts": 0,
        "kc": 0,
    }
)
df.to_csv(output_file_csv, index=False)

param_dicts = df.to_dict(orient="records")
with open(output_file_py, "w") as f:
    f.write("parameters = [\n")
    for row in param_dicts:
        f.write(f"    {row},\n")
    f.write("]\n")
