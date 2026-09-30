"""Experiment: Theorem 7.1 and Observation 7.3 (the multiplier-order threshold).

Claim to test: for EVERY odd prime p, the first depth k at which
P_c(p^2, k) != P_c(p, k) equals T_c(p) = min l * ord_p(lambda) over the
non-critical cycles (infinity if all cycles are critical).
Also test the bound used in the proof of part (3): the number of periodic points
of x^2+c on F_p is at most (p+1)/2 (they lie in the image of f).
"""
import subprocess
import sys

from adt import primes_upto, parse_line, CYCLES_BIN
from exp_primepower import threshold

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
ps = [p for p in primes_upto(PMAX) if p > 2]
agree = mism = blind_ok = blind_bad = big = 0
maxper_ratio = 0.0
for c in (1, 2, 3, 0, -1, -3, 5, -4):
    for i in range(0, len(ps), 60):
        chunk = ps[i:i + 60]
        mods = ",".join(str(p * p) for p in chunk) + "," + ",".join(str(p) for p in chunk)
        out = subprocess.run([CYCLES_BIN, "list", mods, str(c)], capture_output=True, text=True, check=True).stdout
        Z = {n: z for n, cc, z in map(parse_line, out.splitlines())}
        for p in chunk:
            A, B = Z[p], Z[p * p]
            nper = sum(L * m for L, m in A.items())
            maxper_ratio = max(maxper_ratio, nper / ((p + 1) / 2))
            T = threshold(p, c)
            # first k with differing marks P(k) = sum_{d|k} d N_d
            kmax = 1
            for L in list(A) + list(B):
                kmax = max(kmax, L)
            first = None
            for k in range(1, 2 * kmax + 2):
                a = sum(L * m for L, m in A.items() if k % L == 0)
                b = sum(L * m for L, m in B.items() if k % L == 0)
                if a != b:
                    first = k
                    break
            if T is None:
                if first is None:
                    blind_ok += 1
                else:
                    blind_bad += 1
                    print("blind mismatch", c, p, first)
                continue
            if 2 ** min(T, 200) > p:
                big += 1
            if first == T:
                agree += 1
            else:
                mism += 1
                print(f"MISMATCH c={c} p={p} T={T} first={first}")
print(f"odd primes < {PMAX}, 8 probes: agree {agree}, mismatch {mism}; "
      f"of the agreeing/mismatching cases {big} have 2^T > p (so p > 2^T + 1 fails)")
print(f"all-critical primes: {blind_ok} invisible at every depth as predicted, {blind_bad} not")
print(f"max over (p,c) of #periodic points / ((p+1)/2) = {maxper_ratio:.4f}  (must be <= 1)")
