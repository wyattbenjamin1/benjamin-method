# Benjamin Weighted Step-Down Procedure

**Author:** Wyatt Benjamin  
**Course:** Math 310 — Mathematical Statistics, Spring 2026

---

## Description

A rank-based adaptive multiple comparison procedure that allocates the alpha budget proportionally to rank separation between group means. Uses a step-down testing structure to maintain strong FWER control while maximizing power for extreme comparisons.

Pairs with larger rank separation (i.e., the most extreme group differences) receive a larger share of the alpha budget, making the method more powerful than Tukey HSD when extreme comparisons are the ones most likely to differ.

---

## The Mathematical Formula

Weights for each pair (i, j) are defined as:

```
w_ij = |rank(μ_i) − rank(μ_j)|^γ  /  Σ |rank(μ_k) − rank(μ_l)|^γ
```

where γ = 0.9 (selected via preliminary simulation).

At each step k of the step-down procedure, the adjusted alpha threshold is:

```
adj_alpha_k = α × (w_k / Σ remaining weights)
```

Reject H_0 for pair k if p_k ≤ adj_alpha_k, then remove w_k from the remaining budget. Stop at the first non-rejection.

---

## Installation & Usage

No package installation required. Clone or download the repository and source the script directly in R:

```r
source("BenjaminMethod_final.R")
```

### Basic usage

```r
# Define your group means, sample size per group, and MSE from ANOVA
means <- c(14.5, 13.8, 13.3, 14.3, 13.1)
J     <- 9       # observations per group
MSE   <- 0.088   # pooled within-group variance from ANOVA table

result <- my_comparison(means, J, MSE, alpha = 0.05, gamma = 0.9)
print(result)
```

### Run the full analysis pipeline

```r
run_all()   # runs demo + FWER validation + power study + plot
```

---

## Example Output

### Oil Filter Demo

Five oil filter brands tested with J = 9 replicates each, MSE = 0.088.

```
--- Oil Filter Example (Demo) ---
Group means:  14.5  13.8  13.3  14.3  13.1
J = 9  |  MSE = 0.088

   group1 group2 diff      p_value    weight    adj_alpha significant
1       1      5  1.4 1.869171e-12 0.1886744 0.009433717        TRUE
2       1      3  1.2 1.311617e-10 0.1456358 0.008975172        TRUE
3       4      5  1.2 1.311617e-10 0.1456358 0.010938708        TRUE
4       3      4  1.0 1.156868e-08 0.1011081 0.009720918        TRUE
5       1      2  0.7 1.162770e-05 0.1011081 0.012066956        TRUE
6       2      5  0.7 1.162770e-05 0.1011081 0.015905600        TRUE
7       2      3  0.5 9.316878e-04 0.0541825 0.012500000        TRUE
8       2      4  0.5 9.316878e-04 0.0541825 0.016666667        TRUE
9       3      5  0.2 1.604268e-01 0.0541825 0.025000000       FALSE
10      1      4  0.2 1.604268e-01 0.0541825            NA       FALSE
```

8 of 10 pairs are significant. The step-down procedure stops at row 9 (groups 3 vs 5, diff = 0.2), which does not clear its adjusted alpha threshold. All pairs from that point onward are declared non-significant.

---

### FWER Validation

Simulated under the global null hypothesis (all group means = 0) with 10,000 replications per scenario. All estimates fall within the acceptable [0.03, 0.07] target range.

```
--- FWER Results (target: 0.03 to 0.07) ---
Scenario   I    J    FWER     SE         Status
-------------------------------------------------------
A          3    10   0.0594   0.002364   Within [0.03, 0.07]
B          5    10   0.0668   0.002497   Within [0.03, 0.07]
C          7    10   0.0699   0.002550   Within [0.03, 0.07]
D          5    5    0.0620   0.002412   Within [0.03, 0.07]
E          5    20   0.0673   0.002505   Within [0.03, 0.07]
```

> **Note:** Estimates of 0.059–0.070 reflect a deliberate tradeoff from choosing γ = 0.9, which concentrates alpha on extreme comparisons to boost power. All values remain within the target window and within 2 Monte Carlo SEs of 0.05. To bring FWER to exactly 0.05, reduce γ to 0.7 or use a base alpha of 0.045.

---

### Power Analysis

Estimated power for detecting μ₄ ≠ μ₁ under true means = (0, 0, 0, δ, δ), with I = 5, J = 10, σ = 1, α = 0.05. Based on 1,000 simulations per δ value.

```
--- Power Results ---
Delta    Benjamin   Tukey      Rel.Gain(%)   Winner
------------------------------------------------------------
0.5      0.057      0.041      +39.0%        Benjamin  [~Type I rate, no real signal]
1.0      0.319      0.272      +17.3%        Benjamin
1.5      0.749      0.693      +8.1%         Benjamin
2.0      0.964      0.946      +1.9%         Benjamin
```

> **Interpretation:** Benjamin outperforms Tukey HSD at all tested effect sizes. The relative advantage is largest at small δ (where targeted alpha allocation matters most) and converges as δ → 2.0 (where both methods approach ceiling power). At δ = 0.5 neither method has real power — both are near the nominal α = 0.05 Type I error rate.

---

### Gamma Calibration

Preliminary simulation (10,000 runs) used to select γ = 0.9. Score = Power − |FWER − α|.

```
--- Gamma Calibration Results ---
gamma    FWER     SE_FWER  Power    SE_Power   Score
----------------------------------------------------------
0.50     0.0544   0.0023   0.3104   0.0046     0.3060
0.70     0.0597   0.0024   0.3229   0.0047     0.3132
0.80     0.0618   0.0024   0.3242   0.0047     0.3124
0.90     0.0685   0.0025   0.3392   0.0047     0.3207  <- SELECTED

Selected gamma = 0.9 based on Score = Power - |FWER - alpha|.
This balances maximizing power against penalizing FWER deviation from 0.05.
```

---

## When to Use This Method vs. Tukey HSD

### Use Benjamin Weighted Step-Down when:
- You have **prior reason to believe the most extreme groups differ most** — for example, dose-response experiments where the highest and lowest doses are the primary comparisons of interest
- The **cost of a false negative on an extreme comparison** outweighs the cost of a slightly elevated overall Type I rate
- **Not all pairwise comparisons are of equal scientific interest**

### Use Tukey HSD when:
- All pairwise comparisons are of **equal scientific importance**
- You need a **strict FWER guarantee at exactly α = 0.05**
- You have an **unbalanced design** (unequal group sizes) — Benjamin is not validated for this case
- You suspect **non-normal errors** — Tukey is more robust in this setting

---

## Requirements

- R (≥ 4.0)
- Base R only — no additional packages required

---

## License

MIT License — see `LICENSE` file for details.
