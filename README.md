# moment-hoc
This is code for the paper "Optimal Control of Hybrid Systems via Occupation Measures".

## Getting Started
The code requires the SDP solver [MOSEK](https://www.mosek.com/) and, for baseline comparison, the MIQP solver [GUROBI](https://www.gurobi.com/). Academic licenses are available for both of these packages.

The experiment scripts will detect typical installation locations automatically.

## Running the experiments
1. Clone the repository, then open a Julia 1.10 REPL in the cloned directory and type `]` to switch to package mode. Enter `activate .` followed by `instantiate` to download and precompile the dependencies.
2. Run `run_all.sh` in a shell to run all the experiments sequentially. The scripts will output results to the `img` folder.
