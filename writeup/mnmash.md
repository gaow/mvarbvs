---
title: "M&M ASH model and implementation"
date: "August 18, 2016"
output: pdf_document
---
## Notations
| Notation   | Definition    |
|----------|:-------------:|
| $J$ | Number of response variables |
| $N$ | Sample size |
| $Y\mapsto{\rm R}^{N \times J}$ | Matrix of phenotypes |
|$X\mapsto{\rm R}^{N \times P}$| Matrix of genotypes |
| $\beta\mapsto{\rm R}^{P \times J}$ | Matrix of genotypic effect size |
|$\Sigma\mapsto{\rm R}^{J \times J}$ | Covariance matrix of errors |
|$V_p$| Covariance matrix of $\hat{\beta}_p$ |

## Established work
Consider a multivariate, multiple regression problem
\[Y \mid X, \beta, \Sigma \sim N_{N \times J}(X\beta, I_N, \Sigma)\]
The goal is to make inference on $\beta$. 
A Bayesian model for $\beta$ is adopted.
As will be shown in the rest of the document, this framework allows for testing and integrating a wide range of association models to uncover patterns of association with multiple phenotypes.
Consider the case for a single SNP $p$:
\[\hat{\beta}_p \mid \beta_p, V_p \sim N_J(\beta_p, V_p)\]
where $\beta_p$ has a matrix form of `ash` prior:
\[\beta_p \mid \pi, U \sim \sum_{k,l} \pi_{k,l}N_J(0, \omega_lU_k)\]
where $\omega_l$ are given grid values, $U_k$ are given matrices and the mixture components $\pi_k,l$ are learned from data.
This is the `mash` model (Urbut *et al* 2016).

## M&M ASH model
`m&mash` is a multiple regression extension of `mash` via variational inference. 
Without loss of generality we choose $\Sigma = I_J$. This is because our model is 
\begin{eqnarray}
Y &=& X\beta + E \\
E  &\sim & N_{N \times J} (0, I , \Sigma) 
\end{eqnarray}
where we focus on $\beta \Sigma^{-\frac{1}{2}}$, so we can consider the model:
\begin{eqnarray}
Y\Sigma^{-\frac{1}{2}} &=& X \beta \Sigma^{-\frac{1}{2}} + \tilde{E} \\
\tilde{E} &=& E \Sigma^{-\frac{1}{2}} \sim  N_{N \times J} (0, I , I) 
\end{eqnarray}
We assume that $\Sigma$ is known, so we move to a simpler model
\begin{eqnarray}
Y &=& X\beta + E \\
E  &\sim & N_{N \times J} (0, I , I) 
\end{eqnarray}

For the prior of $\beta$, we assume that 
\begin{eqnarray}
\beta_j \sim \sum_t \pi_t N(0,V_t)
\end{eqnarray}
where $$V_t=\omega_k U_l$$ and 
\[\beta=\left\{\begin{array}{c}
\beta_1^T\\
\beta_2^T\\
... \\
\beta_P^T \end{array}\right\} \]
So in this case we consider that $V_t$ with single index $t$.

Our goal is to minimize the K-L divergence:

\begin{eqnarray}
F &=& E_q log\frac{q(\beta)}{p(\beta)p(Y|X,\beta)} \\
  &=& E_q log q(\beta) - E_q log p(\beta) - E_q log p(Y | X, \beta)
\end{eqnarray}

We assume that:

\begin{eqnarray}
q(\beta) = \prod_j q(\beta_j)
\end{eqnarray}
and 
\begin{eqnarray}
q(\beta_j) = \sum_t \alpha_{tj} N(\mu_{tj},S_{tj})
\end{eqnarray}

We denote that 
\begin{eqnarray}
r_j = E_q \beta_j = \sum_t \alpha_{tj}\mu_{tj}
\end{eqnarray}
 and $X^TY = [\phi_1,\cdots, \phi_P]^T$, i.e.
\[ X^TY=\left\{\begin{array}{c}
\phi_1^T\\
\phi_2^T\\
... \\
\phi_P^T \end{array}\right\} \]
\begin{eqnarray} 
E_q logp(Y|X,\beta) &=& -\frac{NJ}{2} log2\pi -\frac{1}{2} E_q \{ tr[(Y-X\beta)^T(Y-X\beta)]\}\\
 &=& -\frac{NJ}{2} log2\pi -\frac{1}{2} E_q \{ tr[Y^TY - Y^TX\beta - \beta^T X^T Y + \beta^T X^T X \beta]\}\\
 &=& c_1 + tr[E_q(\beta^T) X^TY] - \frac{1}{2} E_q[tr(\sum_{i = 1}^P \sum_{j = 1}^P \sum_{k =1}^N 
 \beta_i X_{ki} X_{kj} \beta_j^T) ]\\
 &=& c_1 + tr(\sum_j r_j \phi_j^T) - \frac{1}{2} tr[\sum_{i = 1}^P \sum_{j = 1}^P \sum_{k =1}^N 
  X_{ki} X_{kj}E_q(\beta_i \beta_j^T)] \\
 &= & c_1 + tr(\sum_j r_j \phi_j^T)  - \frac{1}{2} tr(\sum_{i = 1}^P \sum_{j = 1}^P \sum_{k =1}^N 
 X_{ki} X_{kj} r_i r_j^T) \nonumber\\
 &+& \frac{1}{2} tr(\sum_{i = 1}^P \sum_{k =1}^N X_{ki} X_{ki} r_i r_i^T) - \frac{1}{2} tr[\sum_{i = 1}^P
 \sum_{k =1}^N X_{ki} X_{ki} E_q(\beta_i \beta_i^T)]\\
 & = & c_1 + tr(\sum_j r_j \phi_j^T)  - \frac{1}{2} tr(\sum_{i = 1}^P \sum_{j = 1}^P \sum_{k =1}^N 
 X_{ki} X_{kj} r_i r_j^T) + \frac{1}{2} tr(\sum_{i = 1}^P \sum_{k =1}^N X_{ki} X_{ki} r_i r_i^T) \nonumber\\
 &-& \frac{1}{2} tr[\sum_{i = 1}^P \sum_{k =1}^N X_{ki} X_{ki} \sum_t \alpha_{it}(\mu_{it} \mu_{it}^T + S_{it})]
\end{eqnarray}

\begin{eqnarray}
E_q log p(\beta) &=& \sum_j \sum_t \alpha_{tj} \{log\pi_t - \frac{J}{2}log2\pi -\frac{1}{2}log|V_t| -\frac{1}{2}
tr[V_t^{-1} (\mu_{tj} \mu_{tj^T}+ S_{tj})]\}
\end{eqnarray}

\begin{eqnarray}
E_q log q(\beta) = \sum_j \sum_t \alpha_{tj}(log\alpha_{tj} - \frac{J}{2}log2\pi -\frac{1}{2}log|S_{tj}| - 
\frac{J}{2})
\end{eqnarray}

So we can get that
\begin{eqnarray}
\frac{\partial F}{ \partial S_{tj}} &=& \frac{1}{2}\alpha_{tj}(d_jI - S_{tj}^{-1} + V_t^{-1}) \\
\frac{\partial F}{ \partial \mu_{tj}} &=& \alpha_{tj}(-\phi_j +\sum_i\sum_k X_{ki}X_{kj}r_i - d_j
r_j + d_j \mu_{tj} + V_t^{-1}\mu_{tj})\\
\frac{\partial F}{ \partial \alpha_{tj}} &=& - \mu_{tj}^T\phi_j + \sum_i\sum_k X_{ki}X_{kj}\mu_{tj}^Tr_i - d_j
\mu_{tj}^T r_j + \frac{1}{2} d_j tr(\mu_{tj}\mu_{tj}^T + S_{tj}) \nonumber\\
&+& log\frac{\alpha_{tj}}{\pi_t} -\frac{1}{2} log\frac{|S_{tj}|}{|V_t|} + \frac{1}{2} tr[V_t^{-1}(\mu_{tj}\mu_{tj}^T + S_{tj})] - \frac{J}{2} +1
\end{eqnarray}
where $d_j = \sum_k X_{kj}X_{kj}$.

To solve this first order derivative conditions:
\begin{eqnarray}
S_{tj} &=& (d_jI + V_t^{-1})^{-1} \\
\mu_{tj} &=& S_{tj} (\phi_j - \sum_i [X^TX]_{ij} r_i + [X^TX]_{jj}r_j)\\
\alpha_{tj} & \propto & \pi_t \sqrt{\frac{|S_{tj}|}{|V_t|}}\exp \{\frac{1}{2}\mu_{tj}^T S_{tj}^{-1}\mu_{tj}\}\\
\sum_t \alpha_{tj} &=& 1 \nonumber
\end{eqnarray}
or you can also use
\begin{eqnarray}
S_{tj} &=& ([X^TX]_{jj}I + V_t^{-1})^{-1} \\
\mu_{tj} &=& S_{tj} ([Y^TX]_j - \sum_{i \neq j} [X^TX]_{ij} r_i )\\
\alpha_{tj} & \propto & \pi_t \sqrt{\frac{|S_{tj}|}{|V_t|}}\exp \{\frac{1}{2}\mu_{tj}^T S_{tj}^{-1}\mu_{tj}\}\\
\sum_t \alpha_{tj} &=& 1 \nonumber
\end{eqnarray}
So, the equation (24)-(26) or (27)-(29) are the iteration scheme.


