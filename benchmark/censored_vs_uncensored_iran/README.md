# Censored vs uncensored QUAIDS

This benchmark compares QUAIDS with and without the Shonkwiler–Yen censoring
correction on the same Iranian household food-demand data. Both models use
IFGNLS and the built-in expenditure control function (`ivexp`).

| Measure | Censored | Uncensored |
|---|---:|---:|
| Households | 37,656 | 37,656 |
| Goods | 10 | 10 |
| Successful bootstrap replications | 200 | 200 |
| Failed bootstrap replications | 0 | 0 |

All ten own-price (uncompensated) point estimates are negative in both models. Their absolute
values are larger with censoring correction. See the
[results table](results/elasticity_comparison.md) for expenditure and own-price
elasticities with 95% percentile intervals, or the
[final log](results/logs/final_summary.txt). SSBs means sugar-sweetened beverages.

## Data and model

The data come from the authors' Iran added-sugar replication data. Original
group 3 (table sugar) is removed. Budget shares and expenditure are recalculated
for the ten remaining groups. One household has no spending in those groups
and is excluded from both models. Original group numbers are kept in the files.

Both models use the same prices, shares, expenditure, and demographic variables:
`age`, `scale`, `children`, and `region_1` through `region_5`. They use `anot=11`
and no curvature restriction. Prices are prepared in advance and held fixed;
the bootstrap does not include uncertainty from constructing those prices.

The instruments are centered log household income and its square. Survey
weights are used only to center log income. The internal expenditure regression
uses ordinary least squares; demand estimation does not use survey weights.
Its residual enters the demand equations and, with censoring correction, the
participation Probits.

The package bootstrap samples household rows independently with replacement.
It rebuilds the expenditure regression and refits the model in each draw.
Each model uses the same seeds for matching batches, so matching draw IDs refer
to the same household resamples. There is no cluster resampling.

## Files

| File or folder | Contents |
|---|---|
| `data/iran_food_demand_10goods.csv` | Ready-to-use shares, prices, expenditure, demographics, income, and survey weights. |
| `data/group_labels.csv` | Original group numbers and names. |
| `data/dataset_metadata.json` | Data source, preparation steps, columns, and hashes. |
| `scripts/run_bootstrap.py` | Fits one model and saves a batch of draws. |
| `scripts/summarize_bootstrap.py` | Checks and combines batches, then writes tables and file hashes. |
| `results/elasticity_comparison.md` | Point estimates and 95% percentile intervals. |
| `results/elasticity_comparison.csv` | The same results, plus paired own-price differences (censored minus uncensored). |
| `results/*_point_estimates.csv` | All reported full-sample coefficients and elasticities. |
| `results/*_point_elasticities.csv` | Expenditure and own-price elasticities. |
| `results/*_marshallian_elasticities.csv`, `results/*_hicksian_elasticities.csv` | Full price-elasticity matrices; rows are demand groups and columns are price groups. |
| `results/*_bootstrap_draws.csv` | All 200 draws for each model. |
| `results/*_bootstrap_inference.csv` | Bootstrap standard errors and percentile intervals for all reported parameters. |
| `results/raw/`, `results/logs/` | Draws, settings, and original console logs for each batch. |
| `results/bootstrap_manifest.json` | Seeds, successes, failures, workers, and run times. |
| `results/run_environment.json` | Original software versions, source revision, and computer details. |
| `results/checksums.sha256` | SHA-256 hashes, with paths relative to this benchmark folder. |

## Run

Use the `publication` branch containing expenditure control functions for
uncensored models. The saved runs used source revision
`2893f2816aba699e1a7ff5603367f18b0fbc45c6` (package version `1.8.0`);
the original `v1.8.0` release does not include that feature.

From the repository root, install the checked-out package:

```bash
python -m pip install -e .
```

Run both models using five workers. These commands use the saved batch sizes
and seeds: 50 draws with seed `20261019`, then six batches of 25 draws with
seeds `20261020` through `20261025`.

```bash
for MODEL in censored uncensored; do
  python benchmark/censored_vs_uncensored_iran/scripts/run_bootstrap.py \
    --model "$MODEL" --batch first --n-jobs 5

  for FIRST in 51 76 101 126 151 176; do
    SEED=$((20261020 + (FIRST - 51) / 25))
    python benchmark/censored_vs_uncensored_iran/scripts/run_bootstrap.py \
      --model "$MODEL" --batch additional --first-rep "$FIRST" --reps 25 \
      --seed "$SEED" --n-jobs 5
  done
done

python benchmark/censored_vs_uncensored_iran/scripts/summarize_bootstrap.py
```

The models run one after the other. Change `--n-jobs` to set the worker count.
Original worker counts and times are recorded for each batch.

To rebuild the tables from the supplied draws without fitting the models,
run only the summary command. Standard errors use the sample standard deviation
of the draws; intervals use their 2.5th and 97.5th percentiles. Paired difference
intervals use the differences between matching draws.
