---
title: "M&M ASH model and implementation"
date: "August 18, 2016"
output: pdf_document
geometry: margin=0.8in
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
|$V_p\mapsto{\rm R}^{J \times J}$| Covariance matrix of $\hat{\beta}_p$ |
|$U\mapsto{\rm R}^{K \times J \times J}$ | Prior matrices (candidate models) |
|$\omega\mapsto{\rm R}^L$ | grid values for the scales of $U$ |

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
where $\omega_l$ are given grid values, $U_k$ are given matrices and the mixture components $\pi_{k,l}$ are learned from data.
This is the `mash` model (Urbut *et al* 2016).

## M&M ASH model
`m&mash` is a multiple regression extension of `mash` via variational inference. 
Without loss of generality we choose $\Sigma = I_J$. This is because under the assumption 
of independent $Y$, our model is 
\begin{eqnarray}
Y &=& X\beta + E \\
E  &\sim & N_{N \times J} (0, I_N, \Sigma_J) 
\end{eqnarray}
If we scale $\beta$ by $\Sigma^{-\frac{1}{2}}$ we can consider the equivalent model:
\begin{eqnarray}
Y\Sigma^{-\frac{1}{2}} &=& X \beta \Sigma^{-\frac{1}{2}} + E\Sigma^{-\frac{1}{2}} \\
E \Sigma^{-\frac{1}{2}} &\sim&  N_{N \times J} (0, I_N, I_J) 
\end{eqnarray}
Assuming $\Sigma$ is known, here we solve this equivalent model:
\begin{eqnarray}
Y &=& X\beta + E \\
E  &\sim & N_{N \times J} (0, I_N, I_J) 
\end{eqnarray}

For the prior of effect size $\beta$, we assume a multivariate normal prior
(the `mash` prior) on each SNP:
\begin{eqnarray}
\beta_P \sim \sum_t \pi_t N_J(0,W_t)
\end{eqnarray}
where $$W_t=\omega_k U_l$$
Here we consider all $\beta_p$ jointly in our inference,
thus the name `m&m` for multivariate multiple regression:
\[\beta=\left\{\begin{array}{c}
\beta_1^T\\
\beta_2^T\\
... \\
\beta_P^T \end{array}\right\} \]

## Variational inference
The objective function to minimize is the K-L divergence:

\begin{eqnarray}
F &=& E_q log\frac{q(\beta)}{p(\beta)p(Y|X,\beta)} \\
  &=& E_q log q(\beta) - E_q log p(\beta) - E_q log p(Y | X, \beta)
\end{eqnarray}

We take a variational approach assuming that:

\begin{eqnarray}
q(\beta) = \prod_p q(\beta_p)
\end{eqnarray}
where
\begin{eqnarray}
q(\beta_p) = \sum_t \alpha_{pt} N(\mu_{pt},S_{pt})
\end{eqnarray}

We denote that 
\begin{eqnarray}
r_p = E_q \beta_p = \sum_t \alpha_{pt}\mu_{pt}
\end{eqnarray}
 and 
\[ Y^TX=\left\{\begin{array}{c}
\phi_1\\
\phi_2\\
... \\
\phi_P \end{array}\right\} \]
We work out the K-L divergence
\begin{eqnarray} 
E_q logp(Y|X,\beta) &=& -\frac{NJ}{2} log2\pi -\frac{1}{2} E_q \{ tr[(Y-X\beta)^T(Y-X\beta)]\}\\
 &=& -\frac{NJ}{2} log2\pi -\frac{1}{2} E_q \{ tr[Y^TY - Y^TX\beta - \beta^T X^T Y + \beta^T X^T X \beta]\}\\
 &=& c_1 + tr[E_q(\beta^T) X^TY] - \frac{1}{2} E_q[tr(\sum_{i = 1}^P \sum_{p = 1}^P \sum_{k =1}^N 
 \beta_i X_{ki} X_{kp} \beta_p^T) ]\\
 &=& c_1 + tr(\sum_p r_p \phi_p^T) - \frac{1}{2} tr[\sum_{i = 1}^P \sum_{p = 1}^P \sum_{k =1}^N 
  X_{ki} X_{kp}E_q(\beta_i \beta_p^T)] \\
 &= & c_1 + tr(\sum_p r_p \phi_p^T)  - \frac{1}{2} tr(\sum_{i = 1}^P \sum_{p = 1}^P \sum_{k =1}^N 
 X_{ki} X_{kp} r_i r_p^T) \nonumber\\
 &+& \frac{1}{2} tr(\sum_{p = 1}^P \sum_{k =1}^N X_{kp} X_{kp} r_p r_p^T) - \frac{1}{2} tr[\sum_{p = 1}^P
 \sum_{k =1}^N X_{kp} X_{kp} E_q(\beta_p \beta_p^T)]\\
 & = & c_1 + tr(\sum_p r_p \phi_p^T)  - \frac{1}{2} tr(\sum_{i = 1}^P \sum_{p = 1}^P \sum_{k =1}^N 
 X_{ki} X_{kp} r_i r_p^T) + \frac{1}{2} tr(\sum_{p = 1}^P \sum_{k =1}^N X_{kp} X_{kp} r_p r_p^T) \nonumber\\
 &-& \frac{1}{2} tr[\sum_{p = 1}^P \sum_{k =1}^N X_{kp} X_{kp} \sum_t \alpha_{pt}(\mu_{pt} \mu_{pt}^T + S_{pt})]
\end{eqnarray}

\begin{eqnarray}
E_q log p(\beta) &=& \sum_p \sum_t \alpha_{pt} \{log\pi_t - \frac{J}{2}log2\pi -\frac{1}{2}log|W_t| -\frac{1}{2}
tr[W_t^{-1} (\mu_{pt} \mu_{pt}^T+ S_{pt})]\}
\end{eqnarray}

\begin{eqnarray}
E_q log q(\beta) = \sum_p \sum_t \alpha_{pt}(log\alpha_{pt} - \frac{J}{2}log2\pi -\frac{1}{2}log|S_{pt}| - 
\frac{J}{2})
\end{eqnarray}

to optimize:
\begin{eqnarray}
\frac{\partial F}{ \partial S_{pt}} &=& \frac{1}{2}\alpha_{pt}(d_pI - S_{pt}^{-1} + W_t^{-1}) \\
\frac{\partial F}{ \partial \mu_{pt}} &=& \alpha_{pt}(-\phi_p +\sum_{i=1}^P\sum_{k=1}^N X_{ki}X_{kp}r_i - d_p
r_p + d_p \mu_{pt} + W_t^{-1}\mu_{pt})\\
\frac{\partial F}{ \partial \alpha_{pt}} &=& - \mu_{pt}^T\phi_p + \sum_{i=1}^P\sum_{k=1}^N X_{ki}X_{kp}\mu_{pt}^Tr_i - d_p
\mu_{pt}^T r_p + \frac{1}{2} d_p tr(\mu_{pt}\mu_{pt}^T + S_{pt}) \nonumber\\
&+& log\frac{\alpha_{pt}}{\pi_t} -\frac{1}{2} log\frac{|S_{pt}|}{|W_t|} + \frac{1}{2} tr[W_t^{-1}(\mu_{pt}\mu_{pt}^T + S_{pt})] - \frac{J}{2} +1
\end{eqnarray}
where $d_p = \sum_k X_{kp}X_{kp}$.

To solve this first order derivative conditions:
\begin{eqnarray}
S_{pt} &=& (d_pI + W_t^{-1})^{-1} \\
\mu_{pt} &=& S_{pt} (\phi_p - \sum_i [X^TX]_{ip} r_i + [X^TX]_{pp}r_p)\\
\alpha_{pt} & \propto & \pi_t \sqrt{\frac{|S_{pt}|}{|W_t|}}\exp \{\frac{1}{2}\mu_{pt}^T S_{pt}^{-1}\mu_{pt}\}\\
\sum_t \alpha_{pt} &=& 1 \nonumber
\end{eqnarray}
or in another notation:
\begin{eqnarray}
S_{pt} &=& ([X^TX]_{pp}I + W_t^{-1})^{-1} \\
\mu_{pt} &=& S_{pt} ([Y^TX]_p - \sum_{i \neq p} [X^TX]_{ip} r_i )\\
\alpha_{pt} & \propto & \pi_t \sqrt{\frac{|S_{pt}|}{|W_t|}}\exp \{\frac{1}{2}\mu_{pt}^T S_{pt}^{-1}\mu_{pt}\}\\
\sum_t \alpha_{pt} &=& 1 \nonumber
\end{eqnarray}
So, the equation (24)-(26) or (27)-(29) are the iteration scheme. We iterate the procedure until the K-L distance convergences to a constant.
