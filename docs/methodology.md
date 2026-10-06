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

### Endogenous expenditure: integrated `ivexp` control function

When total expenditure is endogenous and excluded instruments $q_h$ are
available, the package estimates the auxiliary OLS equation

$$
\ln m_h = \pi_0 + \pi_p'\ln p_h + \pi_z'z_h + \pi_q'q_h + v_h.
$$

This follows the augmented-regression construction for almost-ideal demand
systems: all exogenous variables in the demand system enter the reduced form,
along with the identifying instruments. For each good, the same generated
residual enters the participation equation with coefficient $\psi_i$,

$$
d_{ih}=\mathbb{I}(X_{ih}'\tau_i + \psi_i v_h + u_{ih}>0),
$$

and the latent demand share with a distinct coefficient $\kappa_i$,

$$
\mu_{ih}=\Phi_{ih}\left(w^Q_{ih}+\kappa_i v_h\right)
+\delta_i\phi_{ih}.
$$

The package reports the reduced-form coefficients and covariance, $R^2$, and
the classical joint F statistic for the excluded instruments. The F statistic
is a relevance diagnostic, not a test of the instruments' exclusion validity.
The instruments must be substantively defensible and excluded from both
structural stages. As with other residual-inclusion estimators, consistency
also depends on a correctly specified triangular reduced form and the
distributional conditions used to justify residual inclusion in the Probit
(Rivers and Vuong 1988). Bootstrapping accounts for the generated regressor;
it cannot repair weak or invalid instruments.

### Optional external control function

With an externally estimated reduced-form residual $v_h$, the package can
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
reduced form inside every bootstrap replication. The internal bootstrap does
this automatically for `ivexp`. It remains disabled for an externally supplied
residual because the package does not know that residual's generating equation.

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

### Analytical inference

The Python option `analytic=True` computes joint sandwich and delta-method
inference after the point estimate. It uses the estimating-equation framework
of [Hardin (2002)](https://www.stata-journal.com/article.html?article=st0018).
Observations must be independent, parameters identified, and the estimating
equations solved to adequate numerical accuracy. `ivexp`, control functions,
and survey-design corrections are not yet supported.

For IFGNLS, the parameter vector $\zeta$ contains active Probit coefficients,
free demand parameters, and the distinct entries of $\Sigma$. Its observation
scores are

$$\psi_t(\zeta)=\begin{pmatrix}s_{\tau,t}\\J_t'\Sigma^{-1}u_t\\\operatorname{vech}(u_tu_t'-\Sigma)\end{pmatrix}.$$

NLS uses identity weights and omits the covariance moments. FGNLS adds the
initial NLS equations and builds the weight moments from their residuals.
For this optional calculation, that initial NLS stage is replayed with the
original settings; the point estimate and solver remain unchanged.

With $\bar\psi=N^{-1}\sum_t\psi_t$, define
$A=N^{-1}\sum_t\partial\psi_t/\partial\zeta'$ and
$B=N^{-1}\sum_t(\psi_t-\bar\psi)(\psi_t-\bar\psi)'$.
The joint covariance is $N^{-1}A^{-1}BA^{-T}$, including cross-Probit and
cross-stage terms. Estimating-equation derivatives include the residual
Hessian terms and are calculated analytically in free coordinates.

For a sample mean $m=N^{-1}\sum_t h_t(\tau)$, its influence is
$\varphi_{m,t}=h_t-m+(\partial m/\partial\tau')\varphi_{\tau,t}$,
where $\varphi_{\zeta,t}=-A^{-1}(\psi_t-\bar\psi)$.
For the existing elasticity function $e(\zeta,m)$,
$\varphi_{e,t}=e_\zeta\varphi_{\zeta,t}+e_m\varphi_{m,t}$ and
$\widehat V_e=N^{-2}\sum_t\varphi_{e,t}\varphi_{e,t}'$.
This includes uncertainty in observed means and in the means of the Probit
CDF, PDF, and linear index. Elasticity derivatives use analytical product,
chain, and quotient rules, checked against independent central differences
of the existing function.
The delta approximation is local: strongly nonlinear elasticities, including
those with a near-zero adjusted share in the denominator, can have S.E.s and
confidence-interval coverage that differ materially from finite-sample results.

`res.analytical` provides the covariance and all three elasticity S.E. arrays.
Its `max_standardized_score` reports the largest absolute mean score divided
by its estimated sampling S.E.; it is a diagnostic, not a stopping rule.

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
7. **Lecocq, S., & Robin, J.-M. (2015)**. Estimating Almost-Ideal Demand Systems with Endogenous Regressors. *The Stata Journal*, 15(2), 554–573.
8. **Rivers, D., & Vuong, Q. H. (1988)**. Limited Information Estimators and Exogeneity Tests for Simultaneous Probit Models. *Journal of Econometrics*, 39(3), 347–366.
