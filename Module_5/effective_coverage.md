# **Effective overage instructions**

## Overview
We have just looked at the statistics regarding the nominal coverage from a sample. The coverage of a sample is important because it reflects the precision of allele frequency estimtes. Precision is similar to variance (it is the inverse of variance!). Larger sample sizes give more precise esimates of allele frequency.

To understand the importance of read depth, consider the analogy of a coin. If we want to estimate the frequency that a coin flip gives a "tail", we would flip that coin multiple times. The mulitple flips of a coin is analagous to the multiple reads at a position in the genome. If you flip the coin ten times, maybe you get 3 heads and 7 tails. That does not mean that the true frequency of tails is 70%! The confidence intervals for our frequency estimate are calculated as: $\sqrt{(pq)/n}$. Therefore, in our coin example, of 3 heads & 7 tails, our best estimate of the probability of the coin is 70% but the 95% CIs are approximately 70% ± 2*SE = 70%±.30% = 40-100%. Not a very precise estimate! If we flipped the coin 1000 times and got 700 heads and 300 tails, then our esimate would be 70%±1.5%, and we would feel much better about our esimates of frequency.

With Pool-Seq, we have to account for two levels of precision. First, we need to account for the precision determined by the read depth. Second, we need to account for the precison determined by the number of flies that went into the pool. The reason that we need to account for both is that we ultimately care about estimating the frequency in the population. The pool of flies is a subsample of the total population (usually), and the reads are a subsample of the flies. When we sample flies, we do so without replacement: a fly can only sampled from once. However, when we sample reads we do so with replacement: a fly contains many copies of its genome throughout its body and therefore reads can arise from these separate copies.

Take an extreme case of a pool of four flies sequenced to 1000X coverage. Clearly, the precision of our allele frequency estimate as it relates to the larger population of flies should be based 4 and not 1000.

Fortunately, there is an easy formula to esimate a statistic that we call "effective coverage." Effective coverage accounts for both types of sampling (flies + reads):

$$\sqrt{(pq)/n}$$


Goals:

Objective:
