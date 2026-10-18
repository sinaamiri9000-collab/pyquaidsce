# pyquaidsce vs quaidsce

This benchmark compares pyquaidsce with Stata's quaidsce_c on the same computer.
Both programs fit the same censored QUAIDS model using IFGNLS, with 2,665
households, 4 goods, 2 demographic variables (`age` and `scale`), and `anot(11)`.
Each program completes 200 bootstrap replications using 8 workers.

| Measure | pyquaidsce | quaidsce_c |
|---|---:|---:|
| Time (seconds) | 122.05 | 10,810.18 |
| Time (minutes) | 2.03 | 180.17 |
| Bootstrap replications | 200 | 200 |

pyquaidsce was about **88.6 times faster** in this run.

## Files

| File | Contents |
|---|---|
| `benchmark_speed_200reps.do` | Stata code that runs both programs and records their times. |
| `benchmark_speed_200reps.log` | Results, bootstrap output, and the timing summary. |
| `food_demand_benchmark_4goods.dta` | Input data used by both programs. |
| `stata_bootstrap_reps_200.dta` | The 200 bootstrap results saved by Stata. |

## Run

Install pyquaidsce, quaidsce_c, and Stata's parallel package. Set Stata's
working directory to this folder, then run:

```stata
do benchmark_speed_200reps.do
```
