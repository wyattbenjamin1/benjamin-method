# Benjamin Weighted Step-Down Procedure

This is a rank-based adaptive multiple comparison procedure for one-way ANOVA.
It allocates the alpha budget proportionally to rank separation between group means,
using a step-down testing structure to maintain strong FWER control while
maximizing power for extreme comparisons.

## Mathematical Formula

For each pair (i, j):
- t-statistic: |x̄_i - x̄_j| / sqrt(2 * MSE / J)
- Weight: w_ij = |rank_i - rank_j|^0.9 / Σ|rank_i - rank_j|^0.9
- Adjusted alpha: α_adj = α * (w_ij / remaining_weight)
- Reject if p-value ≤ α_adj, remove weight, continue step-down

## Usage
```r
source("benjamin_method.R")

means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
J <- 9
MSE <- 0.088
my_comparison(means, J, MSE, alpha = 0.05)
```
## Example Output

FWER Results (all within target range 0.03–0.07):
| Scenario | I | J  | FWER   |
|----------|---|----|--------|
| A        | 3 | 10 | 0.0594 |
| B        | 5 | 10 | 0.0668 |
| C        | 7 | 10 | 0.0699 |
| D        | 5 | 5  | 0.0620 |
| E        | 5 | 20 | 0.0673 |

Power vs Tukey's HSD:
| Delta | Benjamin | Tukey |
|-------|----------|-------|
| 0.5   | 0.057    | 0.041 |
| 1.0   | 0.319    | 0.272 |
| 1.5   | 0.749    | 0.693 |
| 2.0   | 0.964    | 0.946 |

## When to Use This Method vs Tukey's HSD

Use the Benjamin method when:

- You have prior reason to believe extreme group means are most likely to differ
- You want higher power for detecting large mean separations
- Your design is balanced (equal group sizes)
- You are exploring data where a few groups are expected to be outliers
- Power is a priority and a slightly elevated FWER is acceptable 

Use Tukey's HSD when:

- You want a well-established, conservative procedure
- Your design is unbalanced
- You have no prior expectation about which groups differ
- You need results that are easily interpretable to a broad audience
- Strict FWER control at exactly 0.05 is required

## Limitations
- Assumes balanced designs and normality
- Performance under unbalanced groups or heavy-tailed distributions is unknown

## License
MIT
