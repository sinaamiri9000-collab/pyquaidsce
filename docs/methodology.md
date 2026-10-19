# Econometric Methodology & Model Specification

`pyquaidsce` estimates the **Quadratic Almost Ideal Demand System (QUAIDS)** introduced by Banks, Blundell, and Lewbel (1997), incorporating **Ray (1983) demographic scaling** (Poi 2012) and the **Shonkwiler and Yen (1999)** two-step correction for censored consumption.

---

## 1. Demand System Specification

Let $n$ denote the number of goods, $p = (p_1, \dots, p_n)'$ the price vector, $m$ total expenditure, and $z = (z_1, \dots, z_R)'$ a vector of household demographic characteristics.

### Price Indices and Scaling Functions

**Translog price index:**

```math
\ln a(p) = \alpha_0 + \sum_{i=1}^n \alpha_i \ln p_i
+ \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \gamma_{ij} \ln p_i \ln p_j.
```

The user sets $\alpha_0$ through `anot`.

**Cobb-Douglas price aggregator:**

```math
b(p) = \prod_{i=1}^n p_i^{\beta_i}
= \exp\left(\sum_{i=1}^n \beta_i \ln p_i\right).
```

**Ray (1983) demographic scaling factors:**

```math
\bar{m}(z) = 1 + \sum_{r=1}^R \rho_r z_r.
```

```math
c(p, z) = \exp\left(\sum_{i=1}^n \sum_{r=1}^R \eta_{ri} z_r \ln p_i\right).
```

**Real expenditure term:**

```math
D = \ln m - \ln a(p) - \ln \bar{m}(z).
```

### Latent (Uncensored) Budget Share Equations

The latent budget share for good $i$, denoted $w_i^{*}$, is:

```math
\begin{aligned}
w_i^{*} ={}& \alpha_i + \sum_{j=1}^n \gamma_{ij} \ln p_j \\
&+ \left(\beta_i + \sum_{r=1}^R \eta_{ri} z_r\right) D
+ \frac{\lambda_i}{b(p) c(p, z)} D^2.
\end{aligned}
```

---

## 2. Theoretical Restrictions

To be consistent with utility maximization, the following restrictions are imposed during estimation:

1. **Homogeneity of Degree Zero in Prices**:

   ```math
   \sum_{j=1}^n \gamma_{ij} = 0 \quad \forall i
   ```

2. **Slutsky Symmetry**:

   ```math
   \gamma_{ij} = \gamma_{ji} \quad \forall i, j
   ```

3. **Demographic Adding-Up**:

   ```math
   \sum_{i=1}^n \eta_{ri} = 0 \quad \forall r
   ```

4. **Adding-up for Uncensored Models** (`censor=False`):

   ```math
   \sum_{i=1}^n \alpha_i = 1, \quad \sum_{i=1}^n \beta_i = 0, \quad \sum_{i=1}^n \lambda_i = 0
   ```

For censored models, all $n$ equations are estimated. The adding-up restrictions
on $\alpha$, $\beta$, and $\lambda$ above are not imposed, and the transformed
fitted shares are not constrained to sum to one.

---

## 3. Shonkwiler & Yen (1999) Censoring Correction

When some households report zero consumption ($w_i = 0$), the Shonkwiler-Yen approach models participation before fitting the demand equations:

### Step 1: Participation Probits
For each good $i$, estimate a Probit model on the binary participation indicator $d_i = \mathbb{I}(w_i > 0)$:

```math
d_i = \mathbb{I}(X_i' \tau_i + u_i > 0).
```

Here $X_i = [\ln p_1, \dots, \ln p_n, \ln m, z_1, \dots, z_R, 1]'$.

From this, compute the standard normal cumulative distribution $\Phi_i$ and probability density $\phi_i$.

### Step 2: System Estimation on Transformed Shares
The observed budget share equation becomes:

```math
w_i = \Phi_i \cdot w_i^{*} + \delta_i \cdot \phi_i + \varepsilon_i.
```

Here $\delta_i$ is an additional structural parameter to be estimated for each equation.

### Endogenous expenditure: integrated `ivexp` control function

To account for endogenous expenditure, supply instrument columns through
`ivexp`, for example `ivexp=["income"]` in Python or `ivexp(income)` in Stata.
For household $h$, the package first uses ordinary least squares to estimate
log expenditure from log prices, demographics, the instruments $q_h$, and a
constant:

```math
\ln m_h = \pi_0 + \pi_p'\ln p_h + \pi_z'z_h + \pi_q'q_h + v_h.
```

The residual $v_h$ is the difference between observed and predicted log
expenditure. The package adds it to both model stages. Each participation
Probit estimates a coefficient $\psi_i$ on this residual:

```math
d_{ih}=\mathbb{I}(X_{ih}'\tau_i + \psi_i v_h + u_{ih}>0),
```

Each demand equation estimates a separate coefficient $\kappa_i$:

```math
\mu_{ih}=\Phi_{ih}\left(w_{ih}^{*}+\kappa_i v_h\right)
+\delta_i\phi_{ih}.
```

The package reports the expenditure regression coefficients, their covariance
matrix, $R^2$, and a joint F test for the instruments. The F test shows how much
the instruments help explain expenditure after accounting for prices and
demographics. It does not establish that the instruments are valid. Their
choice also needs an economic reason for excluding them from the participation
and demand equations.

The correction assumes that the expenditure equation is correctly specified
and that its residual captures the unobserved link with the two model stages.
For Probit, this also relies on assumptions about the joint distribution of
the expenditure and participation errors (Rivers and Vuong 1988).

With `ivexp`, bootstrap re-estimates the expenditure equation and both model
stages in every replication, so standard errors account for estimating the
residual. Bootstrap does not make a weak or invalid instrument suitable.

### Optional external control function

If you have already estimated an expenditure equation, you can pass its
residual column through `control_function`. The package then adds that residual
to each demand equation:

```math
\mu_{ih}=\Phi_{ih}\left(w_{ih}^{*}+\kappa_i v_h\right)
+\delta_i\phi_{ih}.
```

Each good has its own coefficient $\kappa_i$, reported as `cfcoef_i`. Use
`selection_control_function` to add a residual to the participation Probits
as well. Their residual coefficients are estimated separately from the demand
coefficients.

When calculating price and expenditure elasticities, the package holds the
residual fixed. If your application also allows the residual to change, you
need to calculate that additional effect separately.

For bootstrap inference with an external residual, your own bootstrap code
must re-estimate the expenditure equation and rebuild the residual in every
replication. Use `reps=0` inside each call to this package. Automatic bootstrap
is available with `ivexp`, where the package knows how to estimate the
expenditure equation.

---

## 4. Estimation Methods

The system of $n$ equations is estimated using **Nonlinear Seemingly Unrelated Regression (NLSUR)**:

- **NLS**: Minimizes $\sum_t u_t' u_t$ with identity objective weights.
- **FGNLS**: Calculates $\hat{\Sigma} = \frac{1}{N} \sum_t \hat{u}_t \hat{u}_t'$ from NLS residuals and minimizes $\sum_t u_t' \hat{\Sigma}^{-1} u_t$.
- **IFGNLS**: Iterates the FGNLS estimation and updates $\hat{\Sigma}$ until convergence.

The optimization uses an **analytic Gauss-Newton algorithm** with step-halving (or Levenberg-Marquardt damping), evaluated efficiently via block-diagonal delta transformations.

### Convergence

Let $\theta$ and $\theta^+$ denote the parameters before and after an accepted
inner step, and $Q$ and $Q^+$ their weighted SSR values. With residuals $u_t$
and fitted-value Jacobians $J_t$, $g=\sum_t J_t'\Sigma^{-1}u_t$ and $d_{GN}$
is the undamped Gauss-Newton direction computed at $\theta$.

| Tolerance | Criterion (must be strictly below the tolerance) | Default |
| --- | --- | --- |
| `param_tol` | $C_\theta=\max_j\frac{\lvert\theta_j^+-\theta_j\rvert}{1+\lvert\theta_j\rvert}$ | `1e-5` |
| `objective_tol` | $C_Q=\frac{Q-Q^+}{\max(\lvert Q\rvert,10^{-300})}$ | `1e-7` |
| `gn_tol` | $C_{GN}=\frac{\lvert d_{GN}'g\rvert}{\max(\lvert Q\rvert,10^{-300})}$ | `1e-5` |
| `outer_param_tol` | $R_k=\max_j\frac{\lvert\theta_j^{(k)}-\theta_j^{(k-1)}\rvert}{\lvert\theta_j^{(k-1)}\rvert+10^{-8}}$ | `1e-5` |

After an accepted inner step, any one of the three inner criteria suffices.
If no improving step is accepted, only the GN criterion can certify convergence.
For IFGNLS, $R_k$ compares parameters between outer iterations and must pass
in two consecutive iterations; the final inner solve must also converge.

---

## 5. Demand Elasticities

Elasticities are evaluated at sample means ($\bar{w}, \bar{\ln p}, \bar{\ln m}, \bar{z}$):

### 1. Expenditure (Income) Elasticity:

```math
e_i = \frac{\partial \ln q_i}{\partial \ln m} = 1 + \frac{1}{w_i} \frac{\partial w_i}{\partial \ln m}
```

### 2. Uncompensated (Marshallian) Price Elasticity:

```math
e_{ij}^u = \frac{\partial \ln q_i}{\partial \ln p_j} = -\delta_{ij} + \frac{1}{w_i} \frac{\partial w_i}{\partial \ln p_j}
```

where $\delta_{ij}$ is the Kronecker delta ($\delta_{ij}=1$ if $i=j$, else $0$).

### 3. Compensated (Hicksian) Price Elasticity:
Computed using the Slutsky equation:

```math
e_{ij}^c = e_{ij}^u + e_i \cdot w_j
```

---

## References

1. **Banks, J., Blundell, R., & Lewbel, A. (1997)**. Quadratic Engel Curves and Consumer Demand. *The Review of Economics and Statistics*, 79(4), 527–539.
2. **Caro, J. C., Melo, G., Molina, J. A., & Salgado, J. C. (2021)**. Censored QUAIDS estimation with quaidsce. *Boston College Working Papers in Economics*, 1045.
3. **Deaton, A., & Muellbauer, J. (1980)**. An Almost Ideal Demand System. *The American Economic Review*, 70(3), 312–326.
4. **Poi, B. P. (2012)**. Easy Demand-System Estimation with quaids. *The Stata Journal*, 12(3), 433–446.
5. **Ray, R. (1983)**. Measuring the Costs of Children: An Alternative Approach. *Journal of Public Economics*, 22(1), 89–102.
6. **Shonkwiler, J. S., & Yen, S. T. (1999)**. Two-Step Estimation of a Censored System of Equations. *American Journal of Agricultural Economics*, 81(4), 972–982.
7. **Lecocq, S., & Robin, J.-M. (2015)**. Estimating Almost-Ideal Demand Systems with Endogenous Regressors. *The Stata Journal*, 15(2), 554–573.
8. **Rivers, D., & Vuong, Q. H. (1988)**. Limited Information Estimators and Exogeneity Tests for Simultaneous Probit Models. *Journal of Econometrics*, 39(3), 347–366.
