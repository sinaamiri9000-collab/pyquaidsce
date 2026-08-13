# Controlled Stata–Python runtime comparison

This directory is for a **same-data, same-machine** comparison of `quaidsce`
and `pyquaidsce`. The purpose is to document performance without turning a
single timing observation into a universal speed claim.

## Historical evidence

A July 2026 Stata research log for an 11-good censored QUAIDS model records:

- 37,645 final observations;
- 11 goods;
- 3 demographic variables;
- `method(ifgnls)`;
- log opened: 8 Jul 2026, 19:54:14;
- log closed: 9 Jul 2026, 13:59:49.

The opening-to-closing span is 18 h 05 min 35 s. This is useful historical
wall-clock evidence, but it is **not** treated as the final controlled timing
because the timestamps bracket the whole log rather than only the estimator,
and that run supplied a previously computed `initial(b0)` vector.

The original historical log is retained as
[`historical_result2_stata.log`](historical_result2_stata.log).

The repository's separate `timing_11goods_37k.log` records a ~10 minute Python
run for a problem of similar dimensions, but on a different benchmark dataset.
Those two numbers therefore must not be presented as a direct speed ratio.

## Rules for the public comparison

For a defensible speed ratio, both programs should use:

1. the exact same benchmark data and final rows;
2. the same 11 shares, 11 prices, expenditure variable, and 3 demographics;
3. the same `anot(1.6)` normalization;
4. `ifgnls` in both programs;
5. the same starting convention (the supplied templates use zero/default starts);
6. no bootstrap around the point estimate;
7. the same computer, while reporting CPU, RAM, operating system, Stata version,
   Python version, NumPy/SciPy versions, and BLAS backend;
8. timing only the estimation call, not data import or result export.

When the final benchmark data are available, store a **minimal benchmark file**
containing only the variables required for this model, preferably as a Stata
`.dta` file so both programs read the same stored numbers. Do not publish the
research microdata unless its license and disclosure rules permit redistribution.

## Expected benchmark file

The templates expect:

```text
bench/runtime_comparison/benchmark_11goods_37645.dta
```

with these variables:

- shares: `w1` ... `w11`;
- prices: `tornqvistssb2`, `tornqvistsweetsnack`, `tornqvistsweetmeal`,
  `tornqvisttea`, `tornqvistsoursnack`, `tornqvistfruitveg`,
  `tornqvistcereals`, `tornqvistprotein2`, `tornqvistdairy`,
  `tornqvistoils`, `tornqvistspices`;
- expenditure: `tfexp`;
- demographics: `scale`, `age`, `cfunc`.

The benchmark file should already contain the final 37,645 rows, so filtering
is performed before timing.

## Run Stata

From this directory:

```stata
do stata_runtime.do
```

The script writes `stata_runtime.log` and uses Stata's `timer` command around
only the estimator call.

## Run Python

After installing `pyquaidsce` from the repository root:

```bash
python python_runtime.py
```

The script writes `python_runtime.json` and uses `time.perf_counter()` around
only `quaidsce(...)`.

## Reporting the result

After both runs, report the comparison in this form:

| Program | Data | Method | Start | Wall time | Final N | Log-likelihood |
|---|---|---|---|---:|---:|---:|
| Stata | same benchmark | IFGNLS | zero/default | ... | 37,645 | ... |
| Python | same benchmark | IFGNLS | zero | ... | 37,645 | ... |

Then compute

```text
speedup = Stata wall time / Python wall time
```

The README should describe that number as an **observed speedup in this
controlled benchmark**, not as a universal property of the package.
