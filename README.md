# LRAP
This is a collection of Fortran subroutines for the constrained low-rank approximation of matrices and tensors:
- [low-rank nonnegative matrix and tensor approximation](https://doi.org/10.1007/s40314-023-02211-2) (not to be confused with NMF)
- low-rank nonnegative matrix completion
- low-rank matrix and tensor approximation in the maximum norm

The algorithms are based on the method of [(quasioptimal) alternating projections](https://arxiv.org/abs/2308.16097).

## How to build
Follow the guidelines of the [MARIA-Fortran](https://github.com/sbudzinskiy/MARIA-Fortran) project, which is used for working with low-rank matrices and tensors.

## In the literature
This code was used to carry out the numerical experiments in the following articles:
- Budzinskiy S. Quasioptimal alternating projections and their use in low-rank approximation of matrices and tensors.  
  arXiv: [2308.16097](https://arxiv.org/abs/2308.16097) (2023).
- Budzinskiy S. On the distance to low-rank matrices in the maximum norm.  
  Linear Algebra Appl 688, 44–58. doi:[10.1016/j.laa.2024.02.012](https://doi.org/10.1016/j.laa.2024.02.012) (2024).
- Budzinskiy S. Entrywise tensor-train approximation of large tensors via random embeddings.  
  arXiv: [2403.11768](https://arxiv.org/abs/2403.11768) (2024).
