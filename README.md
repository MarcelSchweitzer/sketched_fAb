# Randomized sketching for Krylov approximations of large-scale matrix functions

## Overview

The code in this repository implements the algorithms proposed in the manuscript [1]. In particular, it allows to reproduce all tables and figures in the manuscript. 

Please note that most scripts are *not optimized with respect to execution time*, as their purpose is comparing and plotting the *error* of different methods, not comparing their run time.

## Reproducing figures and tables from [1]
For reproducing the figures and tables from [1], run the scripts `drive_*.m` at the top level of the folder structure. 

- `drive_convdiff.m`: This script runs the experiment reported in "Section 5.1: Convection-diffusion example" of [1]. It generates the two plots shown in Figure 5.1.
- `drive_network.m`: This script runs the experiment reported in "Section 5.2: Network example" of [1]. It generates the four plots shown in Figure 5.2.
- `drive_qcd.m`: This script runs the first experiment reported in "Section 5.3: Lattice QCD" of [1]. It generates the four plots shown in Figure 5.3.
- `drive_qcd_timings.m`: This script runs the second experiment reported in "Section 5.3: Lattice QCD" of [1]. It generates the LaTeX code for Table 5.1.

## Citation
If you use our code in a scientific publication, we would appreciate your citing:

```bibtex
@article{GuettelSchweitzer2022,
  title={Randomized sketching for {K}rylov approximations of large-scale matrix functions},
  author={G{\"u}ttel, Stefan and Schweitzer, Marcel},
  year={2022},
  url={https://arxiv.org/abs/2208.11447}
}
```

## License


## References
[1] S. Güttel, M. Schweitzer, Randomized sketching for Krylov approximations of large-scale matrix functions, arXiv preprint arXiv:2208.11447 (2022)