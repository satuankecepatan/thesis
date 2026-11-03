# Thesis

hey. so, i have an idea.

we have two things:

- fruit shapes
- the ITS region of ribosomal DNA

do these things correlate with each other? evolutionarily i mean.

yeah, that's my thesis.

---

## the big picture (why this matters for breeding)

the whole framework comes down to **banana breeding**. if we look at how diverse banana fruits are, is that diversity actually connected to the ITS region? and if it *is* correlated, which way does it go?

- **positive correlation ($r > 0$):** bigger ITS differences = bigger fruit differences. basically, the tree of ITS sequences tracks fruit traits. if you want wild new fruit combinations (heterosis), cross parents from distant ITS clades.
- **negative correlation ($r < 0$):** distant ITS lineages end up with similar-looking fruits (convergent selection because humans like specific banana shapes), OR super close clones have wildly different fruits (somatic mutations). for breeding, this means you can pull in disease-resistant genetics from a totally different clade without ruining the fruit type you're aiming for.
- **no correlation ($r \approx 0$):** fruit morphology and ITS are doing their own separate things (thanks to polyploidy like AAA/AAB/ABB, clonal sports, or environment). meaning: don't trust ITS alone to predict what the fruit will look like.

## plan of action

### 1. morphology tree (the fruit side)

- [x] sync up the morphology data (make sure sample IDs match the sequence IDs 1:1)
- [x] analyze the traits (standardize continuous numbers with $z$-scores so big measurements don't drown out small ones, and check PCA / collinearity)
- [x] print out plots until a good dendrogram is found (test distance metrics like Gower or Euclidean + clustering methods like `UPGMA`, `Ward.D2`, and `Complete`)
- [x] **justify the methods statistically:**
  - **cophenetic correlation ($r_{coph}$):** check how well the dendrogram actually represents the raw trait distances (aiming for $r_{coph} \ge 0.80$ to prove the clustering method isn't distorting the data).
  - **bootstrap the clusters (`pvclust`):** run 1,000+ resamplings in R to get Approximately Unbiased (AU) $p$-values ($\ge 95\%$) so we know the fruit clusters aren't just random noise.

### 2. molecular tree (the ITS side)

- [x] redo the whole sequence alignment to IQ-TREE (`MAFFT` with `L-INS-i` for tricky indels $\rightarrow$ `trimAl` to clean up gappy columns $\rightarrow$ `IQ-TREE 2`)
- [x] properly print out the dendrogram (root with an outgroup, plot cleanly with `ggtree` or `dendextend`)
- [x] **justify the methods, again, statistically:**
  - **ModelFinder (`-m MFP`):** let BIC / AICc pick the exact nucleotide substitution model instead of guessing `GTR+G`.
  - **branch support (`-bb 1000 -alrt 1000`):** run 1,000 Ultrafast Bootstraps (UFBoot $\ge 95\%$) + SH-aLRT ($\ge 80\%$) to prove the clades hold up.

### 3. tree vs. tree (the showdown)

- [x] compare the two trees statistically:
  - **tanglegram & entanglement (`dendextend`):** plot both trees face-to-face and measure how tangled the connecting lines are ($0 =$ identical, $1 =$ total chaos).
  - **Baker's Gamma ($\gamma$):** rank correlation of where pairs of bananas branch off in both trees (gives a $-1$ to $+1$ score).
  - **Mantel test (`vegan::mantel`):** correlate the fruit distance matrix directly against the ITS genetic distance matrix (9,999 permutations). **this is the main test that tells us if the correlation is positive ($r_M > 0$), negative ($r_M < 0$), or flat ($r_M \approx 0$).**
  - **Procrustes test (`vegan::protest`) & Robinson-Foulds distance (`phangorn::RF.dist`):** check overall ordination match and see which specific banana cultivars are the "rebels" breaking the correlation.
- [x] interpret the results back into the breeding framework
