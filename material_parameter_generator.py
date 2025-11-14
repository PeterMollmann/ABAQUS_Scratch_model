import pandas as pd
from scipy.stats import qmc
import os
import itertools
import numpy as np


def MaterialParameterGenerator(n_samples=10, sampler_type="sobol", grid_points=None):
    param_ranges = {
        "E": (70e3, 300e3),
        "A": (100, 1500),
        "B": (100, 1700),
        "n": (0.1, 0.8),
        "mu": (0.0, 0.2),
    }

    dim = len(param_ranges)

    # -----------------------------
    # 1) SOBOL
    # -----------------------------
    if sampler_type == "sobol":
        sampler = qmc.Sobol(d=dim, scramble=True, seed=42)
        sample = sampler.random(n=n_samples)

    # -----------------------------
    # 2) LATIN HYPERCUBE
    # -----------------------------
    elif sampler_type == "lhs":
        sampler = qmc.LatinHypercube(d=dim, seed=42)
        sample = sampler.random(n=n_samples)

    # -----------------------------
    # 2) HALTON SEQUENCE
    # -----------------------------
    elif sampler_type == "halton":
        sampler = qmc.Halton(d=dim, seed=42, scramble=True)
        sample = sampler.random(n=n_samples)
    # -----------------------------
    # 2) RANDOM
    # -----------------------------
    elif sampler_type == "random":
        sample = np.random.uniform(low=0, high=1, size=(n_samples, dim))

    # -----------------------------
    # 3) GRID
    # -----------------------------
    elif sampler_type == "grid":
        if grid_points is None:
            raise ValueError(
                "grid_points must be provided as a dict: "
                '{"E": 10, "A": 5, "B": 3, "n": 4, "mu": 2}'
            )

        # Build each dimension's grid
        grids = []
        for key, (low, high) in param_ranges.items():
            if key not in grid_points:
                raise ValueError(f"Missing grid_points entry for parameter '{key}'")
            n_pts = grid_points[key]
            grids.append(np.linspace(low, high, n_pts))

        # Cartesian product
        sample = np.array(list(itertools.product(*grids)))
        print(f"Generated {sample.shape[0]} grid points.")

    else:
        raise ValueError("sampler_type must be 'lhs', 'sobol', or 'grid'.")
    # Scale LHS/Sobol to parameter bounds
    if sampler_type in ["lhs", "sobol", "halton", "random"]:
        l_bounds = [low for low, high in param_ranges.values()]
        u_bounds = [high for low, high in param_ranges.values()]
        sample = qmc.scale(sample, l_bounds, u_bounds)

    # Create DataFrame
    df = pd.DataFrame(sample, columns=param_ranges.keys())
    df.insert(0, "id", [f"{i:05d}" for i in range(1, len(df) + 1)])
    df.insert(1, "rho", 7.8e-9)
    df.insert(3, "nu", 0.3)

    df = df.round(
        {
            "E": 0,
            "A": 0,
            "B": 0,
            "n": 3,
            "mu": 3,
        }
    )

    path = "material_parameters/"
    if not os.path.exists(path):
        os.makedirs(path)
    df.to_csv(path + sampler_type + "_material_parameter_sweep.csv", index=False)
    param_dicts = df.to_dict(orient="records")
    with open(path + sampler_type + "_material_parameter_sweep.py", "w") as f:
        f.write("parameters = [\n")
        for row in param_dicts:
            f.write(f"    {row},\n")
        f.write("]\n")
    return df.to_dict(orient="records")


if __name__ == "__main__":
    # 5*5*5*4*4=2000
    MaterialParameterGenerator(
        sampler_type="grid",
        grid_points={
            "E": 5,
            "A": 5,
            "B": 5,
            "n": 4,
            "mu": 4,
        },
    )
    # n_samples = 2048 to better fit sobol sequence
    MaterialParameterGenerator(n_samples=1024, sampler_type="sobol")
    MaterialParameterGenerator(n_samples=1000, sampler_type="halton")
    MaterialParameterGenerator(n_samples=1024, sampler_type="lhs")
    MaterialParameterGenerator(n_samples=1024, sampler_type="random")
