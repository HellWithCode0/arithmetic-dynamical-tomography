"""
Experiment: when does the fingerprint see a square factor?

Theorem 7.1 of the paper: P_c(p^e, k) = P_c(p, k) whenever f_c^k(x) - x is
separable mod p, i.e. unless some cycle of length l | k has multiplier
lambda = prod 2x_i with lambda^(k/l) = 1 (mod p).  Hence the first depth at
which p^2 differs from p is (for every odd prime p)

    T_c(p) = min over non-critical cycles of  l * ord_p(lambda).

(a) verify the prediction against brute force on Z/p^2 for small p;
(b) distribution of T_c(p): prediction T_c(p) >= p^(1/2 - o(1)) for almost all p.
"""
import os
import subprocess
from statistics import median
from math import log

from adt import primes_upto, mult_order, zset_canon, truncate, CYCLES_BIN, parse_line


def cycles_with_multipliers(p, c):
    f = lambda x: (x * x + c) % p
    state = [0] * p
    out = []
    for s in range(p):
        if state[s]:
            continue
        path, x = {}, s
        while not state[x] and x not in path:
            path[x] = len(path)
            x = f(x)
        if not state[x] and x in path:
            cyc, y = [x], f(x)
            while y != x:
                cyc.append(y)
                y = f(y)
            lam = 1
            for z in cyc:
                lam = lam * 2 * z % p
            out.append((len(cyc), lam))
        for y in path:
            state[y] = 1
    return out


def threshold(p, c):
    T = None
    for L, lam in cycles_with_multipliers(p, c):
        if lam == 0:
            continue
        t = L * mult_order(lam, p)
        T = t if T is None else min(T, t)
    return T


if __name__ == "__main__":
    # (a) brute-force check on p^2
    ps = [p for p in primes_upto(160) if p > 2]
    ok = bad = 0
    for c in (1, 2, 3, 0, -1):
        mods = ",".join(str(p * p) for p in ps) + "," + ",".join(str(p) for p in ps)
        res = subprocess.run([CYCLES_BIN, "list", mods, str(c)], capture_output=True, text=True).stdout
        Z = {n: z for n, cc, z in map(parse_line, res.splitlines())}
        for p in ps:
            T = threshold(p, c)
            if T is None:
                continue
            first = next((k for k in range(1, 4 * p * p) if truncate(Z[p * p], k) != truncate(Z[p], k)), None)
            if first == T:
                ok += 1
            else:
                bad += 1
                print(f"   mismatch c={c} p={p}: predicted T={T}, first differing depth={first}")
    print(f"(a) predicted first depth at which p^2 is visible: {ok} agree, {bad} disagree (p < 160, 5 probes)")

    # (b) distribution
    P = [p for p in primes_upto(20000) if p > 1000]
    print("(b) T_c(p) over primes 1000 < p < 20000")
    for c in (1, 2, 3, 0):
        Ts0 = [threshold(p, c) for p in P]
        blind = [p for p, T in zip(P, Ts0) if T is None]
        Ts = [T for T in Ts0 if T is not None]
        P1 = [p for p, T in zip(P, Ts0) if T is not None]
        ratio = [log(T) / log(p) for T, p in zip(Ts, P1)]
        print(f"  c={c}: primes where every cycle is critical (p^e invisible at ALL depths): {blind}")
        print(f"  c={c}: median log T/log p = {median(ratio):.3f};  frac(T <= 50) = {sum(T <= 50 for T in Ts)/len(P):.4f};"
              f"  frac(T < sqrt p) = {sum(T*T < p for T, p in zip(Ts, P1))/len(P):.4f}")
    Tmin = [min(t for c in (-3, -2, -1, 0, 1, 2, 3) if (t := threshold(p, c)) is not None) for p in P[:600]]
    print(f"  C={{-3..3}} (min over probes), first 600 primes > 1000: median log T/log p = "
          f"{median(log(T)/log(p) for T, p in zip(Tmin, P)):.3f}; min T = {min(Tmin)}")
    # blind primes for small generic c, up to 1e5 via the C tool is too slow here; scan p < 20000
    for c in range(1, 11):
        bl = [p for p in primes_upto(20000)[1:] if threshold(p, c) is None]
        print(f"  all-critical primes p < 20000 for c={c}: {bl}")
