# PSCC-SUC

This repository is a Julia-based research codebase for solving Unit Commitment (UC) problems under uncertainty. It implements multiple approaches (Stochastic Optimization, Robus Optimization, Distributionally Robust Optimization) solved by Benders decomposition, and out-of-sample evaluation tools. 

This repository contains the model implementation, data used for experiments.

# Repository layout

- `src/` — main source code
	- `Optimizer/` — optimization algorithms and cut generators
	- `Struct/` — data structures and parsers
- `Data/` — input datasets for test systems
	- `6bus_JEAS/` and `118_syst_JEAS/` — each contains `generators.csv`, `lines.csv`, `load_distribution_profile.csv`, `maximum_load.csv`
- `resultsPSCC/` — notebooks and results from experiments
- `main.ipynb` — notebook to run the code
- `Project.toml`, `Manifest.toml` — Julia project environment and dependencies