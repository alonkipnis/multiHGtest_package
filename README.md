# multiHGtest

**Testing for differences in survival using multiple hypergeometric P-values**

[![PyPI version](https://badge.fury.io/py/multiHGtest.svg)](https://pypi.org/project/multiHGtest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

The **multiHGtest** package provides nonparametric methods for comparing survival between two populations with sensitivity to *rare and localized* departures from the null hypothesis of equal hazards. Unlike standard tests (e.g., log-rank), this method is designed to detect non-proportional hazard differences that occur at only a few time intervals, whereas the identity of those intervals is unknwon in advance. 

The approach works by:
1. Computing exact hypergeometric P-values at each time interval, testing whether the observed number of events in one group deviates from what is expected under equal hazard rates.
2. Combining these per-interval P-values using **higher criticism** or **Fisher's combination test** to obtain an overall test statistic.

This package implements the methods described in:

> Kipnis, A., Galili, B., and Yakhini, Z. (2026). "Higher criticism for rare and weak non-proportional hazard deviations in survival analysis." *Biometrika*, 113(1).

## Installation

From PyPI:

```bash
pip install multiHGtest
```

From source:

```bash
git clone https://github.com/alonkipnis/multiHGtest_package.git
cd multiHGtest_package
pip install -e .
```

### Requirements

- Python >= 3.8
- numpy >= 1.20
- scipy >= 1.7
- pandas >= 1.3
- [multiple-hypothesis-testing](https://pypi.org/project/multiple-hypothesis-testing/)

## Quick Start

### From survival table data (at-risk counts and events)

```python
import numpy as np
from multiHGtest import hchg_test, fisher_hg_test, hg_test_dashboard

# Survival data: at-risk counts and events at each time interval
Nt1 = np.array([100, 95, 90, 85])  # Group 1 at-risk counts
Nt2 = np.array([100, 92, 88, 82])  # Group 2 at-risk counts
Ot1 = np.array([5, 5, 5, 3])       # Group 1 events
Ot2 = np.array([8, 4, 6, 3])       # Group 2 events

# Higher Criticism test
hc_stat = hchg_test(Nt1, Nt2, Ot1, Ot2, alternative='two-sided')
print(f"HC statistic: {hc_stat:.4f}")

# Fisher combination test
fisher_stat, fisher_pval = fisher_hg_test(Nt1, Nt2, Ot1, Ot2)
print(f"Fisher statistic: {fisher_stat:.4f}, p-value: {fisher_pval:.4f}")

# Full dashboard with all statistics
df, stats = hg_test_dashboard(Nt1, Nt2, Ot1, Ot2)
print(df)
print(stats)
```

### From time-to-event data

If your data is in the standard time-to-event format (e.g., from clinical trials), you can convert it first:

```python
import pandas as pd
from multiHGtest import from_time_to_event_to_survival_table, hchg_test

# Time-to-event data with columns: time, event (1=event, 0=censored), group
data = pd.DataFrame({
    'time':  [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
    'event': [1, 1, 0, 1, 0, 1, 0, 1, 1, 0],
    'group': [1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
})

# Convert to survival table format
result = from_time_to_event_to_survival_table(data)

# Run the HCHG test
hc = hchg_test(
    result['at_risk_1'], result['at_risk_2'],
    result['observed_1'], result['observed_2'],
    alternative='two-sided'
)
print(f"HC statistic: {hc:.4f}")
```

## API Reference

### `hchg_test(Nt1, Nt2, Ot1, Ot2, alternative='two-sided', **kwargs)`

Higher Criticism test of hypergeometric P-values. Returns the HC test statistic.

| Parameter     | Description |
|---------------|-------------|
| `Nt1`, `Nt2`  | At-risk counts per time interval for groups 1 and 2 |
| `Ot1`, `Ot2`  | Event counts per time interval for groups 1 and 2 |
| `alternative` | `'greater'`, `'less'`, or `'two-sided'` (default) |
| `gamma`       | HC threshold parameter (default: 0.4) |
| `randomize`   | Use randomized P-values (default: False) |

### `fisher_hg_test(Nt1, Nt2, Ot1, Ot2, **kwargs)`

Fisher combination test of hypergeometric P-values. Returns `(statistic, p_value)`.

### `hg_test_dashboard(Nt1, Nt2, Ot1, Ot2, **kwargs)`

Comprehensive dashboard returning a DataFrame with per-interval results and a dictionary with HC, Fisher, and MinP statistics.

### `hypergeom_test(k, M, n, N, alternative='greater', randomize=False)`

Exact hypergeometric test for comparing proportions. Supports `'greater'`, `'less'`, and `'two-sided'` alternatives.

### `from_time_to_event_to_survival_table(data, time_col, event_col, group_col)`

Convert time-to-event (survival) data into the at-risk/event count format required by the test functions.

### Legacy aliases

For backward compatibility: `HCHGtest`, `FisherHGtest`, `testHG_dashboard`.

## Running Tests

```bash
pip install -e ".[dev]"
pytest
```

## Citation

If you use this package in your research, please cite:

```bibtex
@article{kipnis2026higher,
  title={Higher criticism for rare and weak non-proportional hazard deviations in survival analysis},
  author={Kipnis, Alon and Galili, Ben and Yakhini, Zohar},
  journal={Biometrika},
  volume={113},
  number={1},
  year={2026}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
