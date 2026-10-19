# Examples

These scripts show the same examples in Python, Stata, and R. Each starts with
a basic model, then adds expenditure endogeneity correction, bootstrap,
parallel workers, and starting values. Comments above each example explain
what to change.

| File | Use |
|---|---|
| [example_quaidsce.py](example_quaidsce.py) | Python examples. |
| [example_quaidsce.do](example_quaidsce.do) | Stata examples using `pyquaidsce`. |
| [example_quaidsce.R](example_quaidsce.R) | R examples using `quaidsce()`. |
| [data/household_demand.csv](data/household_demand.csv) | Shared practice dataset with 800 simulated households. |
| [data/generate_data.py](data/generate_data.py) | Rebuilds the practice dataset with a fixed random seed. |

Install the package and your chosen interface first. See the
[Python installation guide](../docs/getting-started.md#1-installation),
[Stata guide](../stata/README.md), or [R guide](../R/README.md).

From the repository root, run the script for your language:

```bash
python examples/example_quaidsce.py
```

```stata
do examples/example_quaidsce.do
```

```r
source("examples/example_quaidsce.R")
```

The Stata and R scripts download the CSV if it is not found locally. For
Python, keep the script and its `data` folder together. You can also copy one
example at a time into your own work after loading the data and variable lists.

## Data columns

| Columns | Meaning |
|---|---|
| `w1` to `w4` | Spending shares for four goods. Each row sums to one; a zero means no purchase of that good. |
| `p1` to `p4` | Positive prices in the same good order. |
| `expenditure` | Total spending on the four goods. |
| `hh_size` | Household size, from 1 to 5. |
| `urban` | Location indicator: 1 for urban, 0 otherwise. |
| `income` | Simulated household income, used as the instrument in the `ivexp` examples. |

These are simulated data for learning the commands, not survey observations.
The generator creates prices, income, expenditure, and shares, then assigns
some zero purchases and rescales shares to sum to one. The examples use
`anot=8` for this dataset.

Bootstrap examples use 10 replications to keep practice runs short. Increase
the number for research results. To see all available options, read the
[User Guide](../docs/user-guide.md).
