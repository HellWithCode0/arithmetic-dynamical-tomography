# Provenance and verification record

## Newly reproduced on 2026-09-30

The files in `results/reproduced/` were generated from `src/adt/` in a clean invocation of the bundled Python 3.11 runtime. The accompanying `SHA256SUMS` records their byte-level identities.

- `collision_65_119.json`: direct full-cycle enumeration for both moduli and both probes.
- `prime_separation_100000.json`: all odd primes below 100,000 checked using the exact `c=0,-2` formulas.
- `squarefree_100000.json`: all odd squarefree moduli through 100,000 searched using CRT/Burnside synthesis.
- `prime_power_depths.csv`: direct comparisons of `p` and `p^2` for odd primes below 100 and four probes.

The test suite also compares the exact formulas with an independent direct enumerator on ten small primes.

## Supplied historical material

`legacy/code/` and `legacy/results/` are copied from the supplied archive `arithmetic-dynamical-tomography_2026-09-30_REVISION2.zip`. The results include computations reported at substantially larger bounds. They were not regenerated as part of the repository construction and remain labeled legacy for that reason. No paper drafts, referee correspondence, or proposal materials are included.

The source archive reported, among other experiments, searches to `10^6` and `10^7`, prime-power threshold checks on 1,540 `(p,c)` cases, and depth-law data on primes below `10^6`. Readers should consult the relevant legacy script and result file together; these claims are historical provenance statements, not claims of a fresh run.

## Reproduction environment

The clean implementation is pure Python and uses no third-party run-time dependency. Tests use pytest. GitHub Actions runs the tests under Python 3.10 and 3.12 and regenerates the 65/119 output.
