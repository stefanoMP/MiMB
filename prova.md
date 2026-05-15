# Bootstrap-based GLM procedure for gene signatures

This is the outline of the procedure in Pagnotta [2024]. The algorithm produces
a bootstrap matrix $W^*$ (rows are associated with genes) of Wald's test
statistics associated with the coefficients of a GLM

$$Y=\beta_0+\beta_1 C_1+\beta_2 C_2+\cdots+\beta_K C_K$$

where the co-variates are the groups $C_k$, $k=1,2,\ldots,K$, associated
with an *a priori* classification, and $Y$ is the gene expression level.

The biological meaning of groups can be explored through enrichment analysis,
which is fed with the columns of $W^*$ one at a time.

> *After some cleaning, the top elements of each column provide the signature
> associated with the $k^{\mathrm{th}}$ group.*

---

## 1. Generate a synthetic cluster $C_0$

Generate a synthetic cluster $C_0$ with a proportion $1/K$ of samples
from each of the other $K$ clusters.

- Let $n_k$ be the size of cluster $C_k$, for each $k=1,2,\ldots,K$

- The size of $C_0$ is $n_0=\min_k n_k$

- From each $C_k$, sample without replacement $n_0/K$ items to fill $C_0$

- Remove from $C_k$ the elements in $C_k\cap C_0$

---

## 2. Generate a bootstrap clustering $C^*$

Generate a bootstrap clustering $C^*$ from
$C_0,C_1,\ldots,C_K$.

- Let $n_k$ be the size of cluster $C_k$, for each $k=0,1,2,\ldots,K$
- 
  (these sizes are different from those computed in Step 1)

- The size of each bootstrap cluster is

$$m = \frac{2}{3}\min_k n_k$$

- From each $C_k$, sample with replacement $m$ items to fill $C^{*}_{bk}, \qquad k=0,1,\ldots,K$, 

where $b$ indexes the bootstrap iterations.

---

## 3. Mine genes associated with $C_1,C_2,\ldots,C_K$

- Feed the DESeq2 procedure<sup>(1)</sup> with RNA-seq data and the bootstrap clustering

$$C^*_{b0},C^{*}_{b1},\ldots,C^{*}_{bK}$$

- Collect the Wald's test statistics $W_{bkj}$ for each gene
  $j=1,2,\ldots,p$,

  where $p$ is the number of genes, and for each $k^{\mathrm{th}}$ group, with $k=1,2,\ldots,K$.

---

Repeat Steps 1, 2, and 3 for $B$ bootstrap iterations ($b=1,2,\ldots,B$).

Let $W^*$ be the bootstrap matrix<sup>(2)</sup> of Wald's test
statistics defined as

$$
W^{*}_{jk}=\frac{1}{B}\sum_{b=1}^{B} W_{bjk}
$$

where each row corresponds to the $j^{\mathrm{th}}$ gene, and columns are
associated with the clustering

$$
C_1,C_2,\ldots,C_K
$$

---

## Significance filtering

Since the elements of $W^*$ come from Wald's tests (an asymptotically
Gaussian test), a second matrix $P$ can be computed where each entry is

$$
p_{jk}=P\!\left[ Z > W^{*}_{jk} \right]
$$

which represents the $p$-value under the alternative hypothesis

$$
\beta_k>0
$$

Let $\alpha$ be the significance level of the test.

Each entry $W^{*}_{jk}$ is set to $0$ if the corresponding

$$
p_{ij}>\alpha
$$

After this update, $W^*$ becomes a matrix of significant non-vanishing
coefficients.

---

## Notes

1. **DESeq2** [Love et al., 2014] estimates the parameters of a generalized
   linear model to solve an ANOVA-like procedure where the bootstrap
   clustering defines the groups. In this setting,
   $C^{*}_{b0}$ acts as a background (control group).

2. At this stage, the columns of $W^*$ are gene profiles that can be queried
   with an enrichment analysis method to highlight the biology associated
   with each cluster.
