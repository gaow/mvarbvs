---
title: "M&M ASH model and inference"
date: "October 8, 2016"
output: pdf_document
geometry: margin=0.8in
---
## Notations
| Notation   | Definition    |
|----------|:-------------:|
| $J$ | Number of conditions |
| $N$ | Sample size |
| $Y\mapsto{\rm R}^{N \times J}$ | Phenotypes |
|$X\mapsto{\rm R}^{N \times P}$| Genotypes (SNPs) |
| $B\mapsto{\rm R}^{P \times J}$ | Genotypic effect size of $P$ SNPs in $J$ conditions|
| $\beta_p\mapsto{\rm R}^{J}$ | Genotypic effect size of one SNP (a row in $B$) |
|$\Sigma\mapsto{\rm R}^{J \times J}$ | Residual covariance|
|$U\mapsto{\rm R}^{K \times J \times J}$ | Prior matrices (candidate models) |
|$\omega\mapsto{\rm R}^L$ | Grid values, scales of $U$ |

## M&M ASH model
Consider a multivariate, multiple regression problem
\begin{eqnarray*}
Y &=& XB + E \\
E  &\sim & \mathcal{MN} (0, I_N, \Sigma) 
\end{eqnarray*}
Or equivalently,
\begin{eqnarray} \label{eq:model}
Y \mid X, B, \Sigma \sim \mathcal{MN}(XB, I_N, \Sigma)
\end{eqnarray}
The goal is to make inference on effect size $B$. 
We use a unimodal mixture prior (`ash` prior) for $B$ 
\begin{eqnarray} \label{eq:prior}
B \mid \pi, U, \omega \sim \sum_{k}\sum_{l} \pi_{k,l}\mathcal{MN}(0, I_p, \omega_lU_k)
\end{eqnarray}
where $\omega_l$ are given grid values, $U_k$ are given matrices and the mixture components $\pi_{k,l}$ are learned from data.
This is the multivariate + multiple regression extension of the `ash` model (Stephens 2016).

## Variational inference
We use variational inference to solve the model. The objective function to minimize is the K-L divergence:

\begin{eqnarray} \label{eq:kl}
F &=& E_q \log\frac{q(B)}{p(B)p(Y|X,B)} \\
  &=& E_q \log q(B) - E_q \log p(B) - E_q \log p(Y | X, B)
\end{eqnarray}
where

\[B=\left\{\begin{array}{c}
\beta_1^T\\
\beta_2^T\\
... \\
\beta_P^T \end{array}\right\} \]

We take a variational approach assuming that

\begin{eqnarray}
q(B) = \prod_p q(\beta_p)
\end{eqnarray}
where
\begin{eqnarray}
q(\beta_p) = \sum_t \alpha_{pt} N(\mu_{pt},S_{pt})
\end{eqnarray}

We then minimize $F$ and estimate model parameters.

### Known $\Sigma$
For starters we assume $\Sigma$ is known; thus without loss of generality we set $\Sigma = I_J$. 
This is because if we scale $B$ by $\Sigma^{-\frac{1}{2}}$ we can consider the equivalent model:
\begin{eqnarray*}
Y\Sigma^{-\frac{1}{2}} &=& XB\Sigma^{-\frac{1}{2}} + E\Sigma^{-\frac{1}{2}} \\
E \Sigma^{-\frac{1}{2}} &\sim&  \mathcal{MN}(0, I_N, I_J) 
\end{eqnarray*}
and the equivalent model
\begin{eqnarray} \label{eq:model2}
Y \mid X, B, \Sigma \sim \mathcal{MN}(XB, I_N, I_J)
\end{eqnarray}
For (\ref{eq:prior}), let $$V_t=\omega_l U_k$$ we re-parameterize the prior
\begin{eqnarray} \label{eq:prior2}
B \mid \pi, W \sim \sum_t \pi_{t}\mathcal{MN}(0, I_p, V_t)
\end{eqnarray}

We denote that 
\[r_p = E_q \beta_p = \sum_t \alpha_{pt}\mu_{pt} \]
 and 
\[ Y^TX=\left\{\begin{array}{c}
\phi_1\\
\phi_2\\
... \\
\phi_P \end{array}\right\} \]
We work out terms in (\ref{eq:kl})
\begin{eqnarray} 
E_q \log p(Y|X,B) &=& -\frac{NJ}{2} \log2\pi -\frac{1}{2} E_q \{ tr[(Y-XB)^T(Y-XB)]\}\\
 &=& -\frac{NJ}{2} \log2\pi -\frac{1}{2} E_q \{ tr[Y^T Y - Y^T X B - B^T X^T Y + B^T X^T X B]\}\\
 &=& c_1 + tr[E_q(B^T) X^TY] - \frac{1}{2} E_q[tr(\sum_{i = 1}^P \sum_{p = 1}^P \sum_{k =1}^N 
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
E_q \log p(B) &=& \sum_p \sum_t \alpha_{pt} \{\log\pi_t - \frac{J}{2}\log2\pi -\frac{1}{2}\log|V_t| -\frac{1}{2}
tr[V_t^{-1} (\mu_{pt} \mu_{pt}^T+ S_{pt})]\}
\end{eqnarray}

\begin{eqnarray}
E_q \log q(B) = \sum_p \sum_t \alpha_{pt}(\log\alpha_{pt} - \frac{J}{2}\log2\pi -\frac{1}{2}\log|S_{pt}| - 
\frac{J}{2})
\end{eqnarray}

The first derivatives
\begin{eqnarray}
\frac{\partial F}{ \partial S_{pt}} &=& \frac{1}{2}\alpha_{pt}(d_pI - S_{pt}^{-1} + V_t^{-1}) \\
\frac{\partial F}{ \partial \mu_{pt}} &=& \alpha_{pt}(-\phi_p +\sum_{i=1}^P\sum_{k=1}^N X_{ki}X_{kp}r_i - d_p
r_p + d_p \mu_{pt} + V_t^{-1}\mu_{pt})\\
\frac{\partial F}{ \partial \alpha_{pt}} &=& - \mu_{pt}^T\phi_p + \sum_{i=1}^P\sum_{k=1}^N X_{ki}X_{kp}\mu_{pt}^Tr_i - d_p
\mu_{pt}^T r_p + \frac{1}{2} d_p tr(\mu_{pt}\mu_{pt}^T + S_{pt}) \nonumber\\
&+& \log\frac{\alpha_{pt}}{\pi_t} -\frac{1}{2} \log\frac{|S_{pt}|}{|V_t|} + \frac{1}{2} tr[V_t^{-1}(\mu_{pt}\mu_{pt}^T + S_{pt})] - \frac{J}{2} +1
\end{eqnarray}
where $d_p = \sum_k X_{kp}X_{kp}$.
The solutions are
\begin{eqnarray}
S_{pt} &=& (d_pI + V_t^{-1})^{-1} \\
\mu_{pt} &=& S_{pt} (\phi_p - \sum_i [X^TX]_{ip} r_i + [X^TX]_{pp}r_p)\\
\alpha_{pt} & \propto & \pi_t \sqrt{\frac{|S_{pt}|}{|V_t|}}\exp \{\frac{1}{2}\mu_{pt}^T S_{pt}^{-1}\mu_{pt}\}\\
\sum_t \alpha_{pt} &=& 1 \nonumber
\end{eqnarray}
or in another notation
\begin{eqnarray}
S_{pt} &=& ([X^TX]_{pp}I + V_t^{-1})^{-1} \\
\mu_{pt} &=& S_{pt} ([Y^TX]_p - \sum_{i \neq p} [X^TX]_{ip} r_i )\\
\alpha_{pt} & \propto & \pi_t \sqrt{\frac{|S_{pt}|}{|V_t|}}\exp \{\frac{1}{2}\mu_{pt}^T S_{pt}^{-1}\mu_{pt}\}\\
\sum_t \alpha_{pt} &=& 1 \nonumber
\end{eqnarray}
We iterate this procedure until $F$ converges.

### Unknown $\Sigma$
We now treat $\Sigma_{J \times J}$ unknown and estimate it in the VB procedure. When $J$ is large, for example $J = 50$, we will have an additional 2,500 parameters to estimate at each iteration; this is likely problematic. There are two simplifications to this: 1) Assuming $\Sigma$ is diagonal, or 2) Assuming $\Sigma$ is low rank.

#### When $\Sigma$ is Low rank 
Let 
\begin{eqnarray}
\Sigma &=& \sigma^2 I + WW^T 
\end{eqnarray}
where $W$ is low rank. Then removing genetic effect the model equivalent to (\ref{eq:model}) is
\begin{eqnarray} \label{eq:model3}
Z = Y - XB \sim \mathcal{MN}(0, I_N, \sigma^2 I + WW^T )
\end{eqnarray}
Using plug-in estimate of $B$, (\ref{eq:model3}) may be solved in the VB procedure along the lines of Tipping and Bishop (1999) (their EM updates with respect to $W$).
