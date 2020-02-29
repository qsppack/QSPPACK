

# QSP phase factors solvers  

A toolbox for solving phase factors in Quantum signal processing.

## Problems and solvers

Given a real polynomial <img src="http://chart.googleapis.com/chart?cht=tx&chl= \tilde{P}" style="border:none;"> of degree <img src="http://chart.googleapis.com/chart?cht=tx&chl= d" style="border:none;"> with definite parity such that <img src="http://chart.googleapis.com/chart?cht=tx&chl= |\tilde{P}(x)| \le 1, x\in[-1,1]" style="border:none;">, the package contains codes for solving phase factors <img src="http://chart.googleapis.com/chart?cht=tx&chl= \Phi=(\phi_0,\dots,\phi_d)" style="border:none;"> such that
<img src="http://chart.googleapis.com/chart?cht=tx&chl=
\begin{aligned}
        &U_\Phi(x) = e^{\mathrm{i} \phi_0 \sigma_z} \prod_{j=1}^{d} \left[ e^{\mathrm{i} \arccos(x) \sigma_x} e^{\mathrm{i} \phi_j \sigma_z} \right]\\
        &= \left( \begin{array}{cc}
        P(x) & \mathrm{i} Q(x) \sqrt{1 - x^2}\\
        \mathrm{i} Q^*(x) \sqrt{1 - x^2} & P^*(x)
        \end{array} \right),
\end{aligned}
" style="border:none;">
where <img src="http://chart.googleapis.com/chart?cht=tx&chl= \Re[P]=\tilde{P}" style="border:none;">.

The package contains two kinds of solvers:

- Optimization-based solver 
- Direct solver (namely the GSLW method and the Haah method)

The package also contains an implementation of the Remez algorithm for finding polynomial approximation.

Applications have been solved by these solvers:

- Hamiltonian simulation
- Eigenstate filter
- Matrix inversion

## References

- [Y. Dong, X. Meng, K. B. Whaley, and L. Lin. Efficient Phase Factor Evaluation in Quantum Signal Processing. arXiv: 2002.11649](https://arxiv.org/abs/2002.11649)
- [A. Gilyén, Y. Su, G. H. Low, and N. Wiebe. Quantum singular value transformation and beyond: exponential improvements for quantum matrix arithmetics. In Proceedings of the 51st Annual ACM SIGACT Symposium on Theory of Computing, pages 193–204, 2019](https://dl.acm.org/doi/10.1145/3313276.3316366)
- [J. Haah. Product decomposition of periodic functions in quantum signal processing.Quantum, 3:190, 2019](https://quantum-journal.org/papers/q-2019-10-07-190/)

  
## Citing our work
If you find our work useful or you use our work in your own project, please consider to cite our work.

## The Authors

We hope that the package is useful for your application. If you have any bug reports or comments, please feel free to email one of the software authors:

* Xiang Meng, mengxianglgal@gmail.com

* Yulong Dong, dongyl@berkeley.edu

  

## Installation

`>> startup`

`>> cd Examples`

`>> test_HS`




