from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .core import canon, periodic_cycles, primes_below, solvable_fingerprint


def squarefree_factors(n: int):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return None
            out.append(d)
        d += 1
    if n > 1:
        out.append(n)
    return out


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def collision(args):
    a = {c: canon(periodic_cycles(65, c)) for c in (0, -2)}
    b = {c: canon(periodic_cycles(119, c)) for c in (0, -2)}
    payload = {"claim": "65 and 119 collide for probes {0,-2}", "verified": a == b,
               "65": {str(k): v for k, v in a.items()}, "119": {str(k): v for k, v in b.items()}}
    write_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True))


def separation(args):
    seen, collisions = {}, []
    for p in primes_below(args.bound):
        if p == 2:
            continue
        fp = solvable_fingerprint((p,))
        if fp in seen:
            collisions.append([seen[fp], p])
        else:
            seen[fp] = p
    payload = {"bound_exclusive": args.bound, "odd_primes_checked": len(seen) + len(collisions),
               "collision_count": len(collisions), "collisions": collisions,
               "scope": "finite verification supporting the unconditional prime-separation theorem"}
    write_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True))


def squarefree(args):
    seen, collisions, checked = {}, [], 0
    for n in range(3, args.bound + 1, 2):
        fs = squarefree_factors(n)
        if fs is None:
            continue
        checked += 1
        fp = solvable_fingerprint(tuple(fs))
        if fp in seen:
            collisions.append([seen[fp], n])
        else:
            seen[fp] = n
    payload = {"bound_inclusive": args.bound, "odd_squarefree_checked": checked,
               "collision_count": len(collisions), "collisions": collisions,
               "scope": "new finite rerun; compare legacy/results for the supplied larger-bound search"}
    write_json(args.output, payload)
    print(json.dumps(payload, sort_keys=True))


def depth(args):
    rows = []
    for p in primes_below(args.prime_bound):
        if p == 2:
            continue
        for c in args.probes:
            zp, zp2 = periodic_cycles(p, c), periodic_cycles(p * p, c)
            max_depth = max([1, *zp, *zp2]) * 2
            first = next((k for k in range(1, max_depth + 1)
                          if sum(d*m for d,m in zp.items() if k%d == 0) !=
                             sum(d*m for d,m in zp2.items() if k%d == 0)), None)
            rows.append({"p": p, "c": c, "first_visible_depth": first})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["p", "c", "first_visible_depth"])
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"rows": len(rows), "output": str(args.output)}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(prog="adt-reproduce")
    sub = parser.add_subparsers(required=True)
    p = sub.add_parser("collision-65-119"); p.add_argument("--output", type=Path, required=True); p.set_defaults(fn=collision)
    p = sub.add_parser("prime-separation"); p.add_argument("--bound", type=int, default=100000); p.add_argument("--output", type=Path, required=True); p.set_defaults(fn=separation)
    p = sub.add_parser("squarefree-search"); p.add_argument("--bound", type=int, default=100000); p.add_argument("--output", type=Path, required=True); p.set_defaults(fn=squarefree)
    p = sub.add_parser("prime-power-depths"); p.add_argument("--prime-bound", type=int, default=100); p.add_argument("--probes", type=int, nargs="+", default=[0, -2, 1, 2]); p.add_argument("--output", type=Path, required=True); p.set_defaults(fn=depth)
    args = parser.parse_args(); args.fn(args)


if __name__ == "__main__":
    main()
