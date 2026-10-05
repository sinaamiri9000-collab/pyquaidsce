"""Rebuild one natural failed bootstrap draw with the original seed.

    python -m validation.reproduce --failed-draw boot-large-base-16001:30

This does not reroll the seed or change starts, solvers, or stopping rules.
"""
import argparse
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from .run import ROOT, fixture, write_json, write_csv
from .cases import BASE, scenarios
from pyquaidsce import quaidsce


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--failed-draw", required=True, help="batch ID:original replication number")
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "results" / "reproductions")
    args = parser.parse_args()
    try:
        id, replication = args.failed_draw.rsplit(":", 1)
        replication = int(replication)
        case = next(c for c in scenarios(100) if c.id == id)
    except (ValueError, StopIteration):
        parser.error("unknown batch ID or replication number")
    count = case.kwargs.get("reps", 0)
    seed = case.kwargs.get("seed")
    if not 1 <= replication <= count or seed is None:
        parser.error("replication must belong to a seeded bootstrap batch")
    original = pd.read_stata(ROOT / "bench" / "small4.dta")
    df, changes = fixture(original, case.fixture)
    kw = dict(BASE, **changes)
    kw.update(case.kwargs)
    child_seed = int(np.random.SeedSequence(seed).generate_state(count, dtype=np.uint32)[replication-1])
    indices = np.random.default_rng(child_seed).integers(0, len(df), size=len(df))
    sample = df.iloc[indices].reset_index(drop=True)
    kw.update(reps=0, verbose=False, sigma_tol=kw.get("boot_sigma_tol", 1e-5))
    started = time.perf_counter()
    with threadpool_limits(limits=1, user_api="blas"):
        res = quaidsce(sample, **kw)
    row = dict(batch=id, original_replication=replication, parent_seed=seed, child_seed=child_seed,
               kwargs=kw, converged=res.converged, n_gn=res.n_gn, n_outer=res.n_outer,
               llf=res.llf, notes=res.notes, elapsed_seconds=time.perf_counter()-started,
               sampled_original_row_indices=indices, theta=res.theta, sigma=res.sigma)
    args.output.mkdir(parents=True, exist_ok=True)
    write_json(args.output / f"{id}-rep{replication}.json", row)
    write_csv(args.output / f"{id}-rep{replication}.csv",
              [dict(name=name, estimate=b, se=res.se[i]) for i, (name, b) in enumerate(zip(res.names, res.b))])
    print(f"{id}, original draw {replication}: converged={res.converged}, n_outer={res.n_outer}")
    return int(not res.converged)


if __name__ == "__main__":
    sys.exit(main())
