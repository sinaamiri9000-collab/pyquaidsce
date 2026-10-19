# Censored vs uncensored QUAIDS

This benchmark shows how correcting for zero expenditure affects estimated
food-demand elasticities. It fits QUAIDS with and without the Shonkwiler–Yen
censoring correction to the same data: 37,656 households in Iran and ten food
groups.

Both models use IFGNLS, demographic controls, and an expenditure control function
(`ivexp`). The instruments are log household income and its square. Each model
uses 200 bootstrap replications, with households sampled independently with
replacement. The expenditure regression is estimated again in every replication.
Prices are held fixed during the bootstrap.

All ten own-price point estimates are negative in both models. They are larger
in absolute value with the censoring correction. The log reports expenditure
and own-price elasticities with 95% bootstrap confidence intervals.

## Files

| File | Contents |
|---|---|
| `iran_food_demand.csv` | The dataset used by both models. |
| `benchmark.py` | Runs both models and writes the results. |
| `benchmark.log` | Model settings and elasticity results. |
| `bootstrap_results.zip` | All 200 bootstrap draws and full-sample estimates for each model, plus the estimation settings. |

## Data

`w1`–`w10` are budget shares and `p1`–`p10` are the corresponding prices.
`expenditure` is total spending on these ten groups. Food group names appear
in the log and at the beginning of the script, in the same order as the columns.

The demographic variables are `age`, `scale` (household size), `children`, and
five regional indicators (`region_1`–`region_5`). The dataset also includes
household `income` and `sample_weight`. Survey weights are used to center log
income; the model equations are estimated without survey weights.

## Run

Install the package from the repository root:

```bash
python -m pip install -e .
```

Then run:

```bash
python benchmark/censored_vs_uncensored_iran/benchmark.py
```

This estimates both models using five workers and updates `benchmark.log` and
`bootstrap_results.zip`. Change `WORKERS` at the top of the script to use a
different number of workers.

To rebuild only the log from the saved results, without estimating the models:

```bash
python benchmark/censored_vs_uncensored_iran/benchmark.py --summarize
```
