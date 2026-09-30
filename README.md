# Arithmetic Dynamical Tomography: reproducibility repository

This repository accompanies **“Arithmetic Dynamical Tomography: Reconstruction Obstructions and Prime-Power Detection.”** It contains a small, dependency-free Python implementation of the central finite computations, tests, machine-readable reproduced outputs, and an archival copy of the supplied earlier computational materials.

For `f_c(x)=x^2+c` on `Z/nZ`, the code records the periodic functional graph as a finite `Z`-set: a mapping from cycle length to number of cycles. It implements the Burnside product used by the Chinese remainder theorem and the closed formulas for the exactly solvable probes `c=0,-2`.

## Quick start

Python 3.10 or later is required. No run-time packages are needed.

```sh
python -m pip install -e ".[test]"
pytest -q
adt-reproduce collision-65-119 --output results/reproduced/collision_65_119.json
adt-reproduce prime-separation --bound 100000 --output results/reproduced/prime_separation_100000.json
adt-reproduce squarefree-search --bound 100000 --output results/reproduced/squarefree_100000.json
adt-reproduce prime-power-depths --prime-bound 100 --probes -3 -2 -1 0 1 2 3 --output results/reproduced/prime_power_depths.csv
```

All commands are deterministic. JSON output uses sorted keys; CSV rows have a fixed iteration order.

## What each command establishes

| command | computation | status of the resulting evidence |
|---|---|---|
| `collision-65-119` | Directly enumerates cycles modulo 65 and 119 for `c=0,-2` | Exact verification of the stated collision |
| `prime-separation` | Uses the proved closed formulas and checks for duplicate fingerprints among odd primes below a chosen bound | Finite regression check; the theorem itself is proved in the article |
| `squarefree-search` | Searches all odd squarefree integers through the chosen bound using the closed formulas and CRT/Burnside multiplication | Newly reproduced finite search |
| `prime-power-depths` | Computes `min ell*ord_p(lambda)` from cycles modulo `p` and independently enumerates the complete cycle inventory modulo `p^2` | Newly reproduced exact small-bound threshold check |

The checked-in `results/reproduced/` files were generated from the clean implementation in `src/`. Larger searches supplied with the project, including bounds up to `10^6` and `10^7`, are in `legacy/results/`. They are retained as historical supplied data and are not represented as fresh reruns by this repository revision.

## Repository layout

- `src/adt/`: portable computational core and command-line interface.
- `tests/`: independent direct-enumeration and algebraic regression tests.
- `results/reproduced/`: outputs generated during this repository build.
- `legacy/code/`: complete supplied computational scripts, including the C enumerator and figure scripts.
- `legacy/results/`: complete supplied result files, unchanged.
- `PROVENANCE.md`: precise distinction between supplied and newly generated artifacts.

The legacy C program can be compiled on a POSIX system with a C compiler; consult `legacy/code/gen_data.sh`. Those scripts are included for auditability and extension of the searches. The clean Python commands above are the supported portable interface.

## Scope and interpretation

Computations test instances and expose exact examples; they do not replace the manuscript’s proofs. In particular, the absence of a collision below a finite bound is not an injectivity theorem. The large numerical files in `legacy/results/` were supplied with the earlier project and should be cited with their bounds and provenance intact.

No license has been assigned because no licensing instruction accompanied the source materials. Citation metadata names Aryaman Katoch as the author.

For prime-power output, an empty threshold has one precise meaning: every cycle modulo `p` is critical and the complete cycle inventories modulo `p` and `p^2` are equal. Thus all marks agree at every depth. The implementation does not infer infinity merely because a finite search cutoff was reached.
