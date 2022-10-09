# optimal_tiebreaker_designs

If you use this code, please cite our [arXiv paper](https://arxiv.org/abs/2202.12511) (Li and Owen, 2022).

Currently we have only implemented the fixed-$x$ optimal monotone designs from Li and Owen (2022),
which we believe corresponds to the setting of greatest interest to practitioners.

All the code is currently in a single R source file, `optimal_fixed_x.R`. 
The main functions are `get_p_dag()`, which computes the essentially two-level designs $p_{\max}^{&dagger;}$ and $p_{\min}^{&dagger;}$ from the text,
and `get_opt_2_level()`, which computes an essentially optimal two-level design, given the treatment fraction and short-term gain constraints.
The efficiency function used is currently the scaled D-optimality criterion `Eff` defined in our paper.
It corresponds to the scaled asymptotic variance of $\hat{\beta}_3$ in the two-line regression model referenced in the paper.
You can define your own (inverse) efficiency criterion, following the definition of `inv_eff()` at the top of `optimal_fixed_x.R`. 
The Head Start example of Section 5.3 is provided in full. The data we used is in `head_start.csv` and is courtesy of [Cattaneo et al. (2017)](https://github.com/rdpackages-replication/CTV_2017_JPAM).



