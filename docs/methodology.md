# Econometric Methodology & Model Specification

`pyquaidsce` estimates the **Quadratic Almost Ideal Demand System (QUAIDS)** introduced by Banks, Blundell, and Lewbel (1997), incorporating **Ray (1983) demographic scaling** (Poi 2012) and the **Shonkwiler and Yen (1999)** two-step correction for censored consumption.

---

## 1. Demand System Specification

Let $n$ denote the number of goods, $p = (p_1, \dots, p_n)'$ the price vector, $m$ total expenditure, and $z = (z_1, \dots, z_R)'$ a vector of household demographic characteristics.

### Price Indices and Scaling Functions

1. **Translog Price Index $\ln a(p)$**:
   $$\ln a(p) = \alpha_0 + \sum_{i=1}^n \alpha_i \ln p_i + \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \gamma_{ij} \ln p_i \ln p_j$$
   where $\alpha_0$ is set by the user (`anot`).

2. **Cobb-Douglas Price Aggregator $b(p)$**:
   $$b(p) = \prod_{i=1}^n p_i^{\beta_i} = \exp\left(\sum_{i=1}^n \beta_i \ln p_i\right)$$

3. **Ray (1983) Demographic Scaling Factors**:
   $$\bar{m}(z) = 1 + \sum_{r=1}^R \rho_r z_r$$
   $$c(p, z) = \exp\left(\sum_{i=1}^n \sum_{r=1}^R \eta_{ri} z_r \ln p_i\right)$$

4. **Deflated Real Expenditure Term $D$**:
   $$D = \ln m - \ln a(p) - \ln \bar{m}(z)$$

### Latent (Uncensored) Budget Share Equations

The latent budget share for good $i$, denoted $w_i^*$, is given by:
$$w_i^* = \alpha_i + \sum_{j=1}^n \gamma_{ij} \ln p_j + \left(\beta_i + \sum_{r=1}^R \eta_{ri} z_r\right) D + \frac{\lambda_i}{b(p) c(p, z)} D^2$$

---

## 2. Theoretical Restrictions

To be consistent with utility maximization, the following restrictions are imposed during estimation:

1. **Homogeneity of Degree Zero in Prices**:
   $$\sum_{j=1}^n \gamma_{ij} = 0 \quad \forall i$$

2. **Slutsky Symmetry**:
   $$\gamma_{ij} = \gamma_{ji} \quad \forall i, j$$

3. **Demographic Adding-Up**:
   $$\sum_{i=1}^n \eta_{ri} = 0 \quad \forall r$$

4. **Adding-up for Uncensored Models** (`censor=False`):
   $$\sum_{i=1}^n \alpha_i = 1, \quad \sum_{i=1}^n \beta_i = 0, \quad \sum_{i=1}^n \lambda_i = 0$$

*(In censored models, adding-up does not hold on observed shares because of the Shonkwiler–Yen transformation, so all $n$ equations are estimated).*

---

## 3. Shonkwiler & Yen (1999) Censoring Correction

When households report zero consumption for some goods ($w_i = 0$), estimating $w_i^*$ directly leads to selection bias. The Shonkwiler–Yen two-step approach addresses this:

### Step 1: Participation Probits
For each good $i$, estimate a Probit model on the binary participation indicator $d_i = \mathbb{I}(w_i > 0)$:
$$d_i = \mathbb{I}(X_i' \tau_i + v_i > 0)$$
where $X_i = [\ln p_1, \dots, \ln p_n, \ln m, z_1, \dots, z_R, 1]'$.

From this, compute the standard normal cumulative distribution $\Phi_i$ and probability density $\phi_i$.

### Step 2: System Estimation on Transformed Shares
The observed budget share equation becomes:
$$w_i = \Phi_i \cdot w_i^* + \delta_i \cdot \phi_i + \varepsilon_i$$
where $\delta_i$ is an additional structural parameter to be estimated for each equation.

### Optional external control function

With an externally estimated reduced-form residual $v_h$, version 1.2.0 can
augment the latent demand share as

$$
\mu_{ih}=\Phi_{ih}\left(w^Q_{ih}+\kappa_i v_h\right)
+\delta_i\phi_{ih}.
$$

Each equation has its own unrestricted $\kappa_i$ (`cfcoef_i`). A residual may
also enter the participation Probits through
`selection_control_function`; those Probit coefficients are distinct from
$\kappa_i$. For price or expenditure perturbations the residual is held fixed,
so the package reports a structural derivative conditional on the supplied
residual. If the reduced form itself changes under the perturbation, its
additional $dv/dt$ term must be handled by the empirical design outside the
generic package.

Because the residual is generated, valid final inference must re-estimate its
reduced form inside every bootstrap replication. The internal bootstrap is
therefore disabled whenever this extension is active.

---

## 4. Estimation Methods

The system of $n$ equations is estimated using **Nonlinear Seemingly Unrelated Regression (NLSUR)**:

- **NLS**: Minimizes $\sum_t u_t' u_t$ assuming an identity error covariance ($\Sigma = I$).
- **FGNLS**: Calculates $\hat{\Sigma} = \frac{1}{N} \sum_t \hat{u}_t \hat{u}_t'$ from NLS residuals and minimizes $\sum_t u_t' \hat{\Sigma}^{-1} u_t$.
- **IFGNLS**: Iterates the FGNLS estimation and updates $\hat{\Sigma}$ until both parameter estimates and covariance matrix converge.

The optimization uses an **analytic Gauss-Newton algorithm** with step-halving (or Levenberg-Marquardt damping), evaluated efficiently via block-diagonal delta transformations.

---

## 5. Demand Elasticities

Elasticities are evaluated at sample means ($\bar{w}, \bar{\ln p}, \bar{\ln m}, \bar{z}$):

### 1. Expenditure (Income) Elasticity:
$$e_i = \frac{\partial \ln q_i}{\partial \ln m} = 1 + \frac{1}{w_i} \frac{\partial w_i}{\partial \ln m}$$

### 2. Uncompensated (Marshallian) Price Elasticity:
$$e_{ij}^u = \frac{\partial \ln q_i}{\partial \ln p_j} = -\delta_{ij} + \frac{1}{w_i} \frac{\partial w_i}{\partial \ln p_j}$$
where $\delta_{ij}$ is the Kronecker delta ($\delta_{ij}=1$ if $i=j$, else $0$).

### 3. Compensated (Hicksian) Price Elasticity:
Computed using the Slutsky equation:
$$e_{ij}^c = e_{ij}^u + e_i \cdot w_j$$

---

## References

1. **Banks, J., Blundell, R., & Lewbel, A. (1997)**. Quadratic Engel Curves and Consumer Demand. *The Review of Economics and Statistics*, 79(4), 527–539.
2. **Caro, J. C., Melo, G., Molina, J. A., & Salgado, J. C. (2021)**. Censored QUAIDS estimation with quaidsce. *Boston College Working Papers in Economics*, 1045.
3. **Deaton, A., & Muellbauer, J. (1980)**. An Almost Ideal Demand System. *The American Economic Review*, 70(3), 312–326.
4. **Poi, B. P. (2012)**. Easy Demand-System Estimation with quaids. *The Stata Journal*, 12(3), 433–446.
5. **Ray, R. (1983)**. Measuring the Costs of Children: An Alternative Approach. *Journal of Public Economics*, 22(1), 89–102.
6. **Shonkwiler, J. S., & Yen, S. T. (1999)**. Two-Step Estimation of a Censored System of Equations. *American Journal of Agricultural Economics*, 81(4), 972–982.
