

# QSP phase factors solvers  

A toolbox for solving phase factors in quantum signal processing.

## Official Website
https://qsppack.gitbook.io/qsppack/
You may find useful tutorials and examples in this website.

## Problems and solvers

Given a real polynomial <img src="http://chart.googleapis.com/chart?cht=tx&chl= f" style="border:none;"> of degree <img src="http://chart.googleapis.com/chart?cht=tx&chl= d" style="border:none;"> with definite parity such that <img src="http://chart.googleapis.com/chart?cht=tx&chl= |f(x)| \le 1, x\in[-1,1]" style="border:none;">, the package contains codes for solving phase factors <img src="http://chart.googleapis.com/chart?cht=tx&chl= \Phi=(\phi_0,\dots,\phi_d)" style="border:none;"> such that

<a href="https://www.codecogs.com/eqnedit.php?latex=U_\Phi(x)&space;=&space;e^{\mathrm{i}&space;\phi_0&space;\sigma_z}&space;\prod_{j=1}^{d}&space;\left[&space;e^{\mathrm{i}&space;\arccos(x)&space;\sigma_x}&space;e^{\mathrm{i}&space;\phi_j&space;\sigma_z}&space;\right]&space;=&space;\left(&space;\begin{array}{cc}&space;P(x)&space;&&space;\mathrm{i}&space;Q(x)&space;\sqrt{1&space;-&space;x^2}\\&space;\mathrm{i}&space;Q^*(x)&space;\sqrt{1&space;-&space;x^2}&space;&&space;P^*(x)&space;\end{array}&space;\right)," target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_\Phi(x)&space;=&space;e^{\mathrm{i}&space;\phi_0&space;\sigma_z}&space;\prod_{j=1}^{d}&space;\left[&space;e^{\mathrm{i}&space;\arccos(x)&space;\sigma_x}&space;e^{\mathrm{i}&space;\phi_j&space;\sigma_z}&space;\right]&space;=&space;\left(&space;\begin{array}{cc}&space;P(x)&space;&&space;\mathrm{i}&space;Q(x)&space;\sqrt{1&space;-&space;x^2}\\&space;\mathrm{i}&space;Q^*(x)&space;\sqrt{1&space;-&space;x^2}&space;&&space;P^*(x)&space;\end{array}&space;\right)," title="U_\Phi(x) = e^{\mathrm{i} \phi_0 \sigma_z} \prod_{j=1}^{d} \left[ e^{\mathrm{i} \arccos(x) \sigma_x} e^{\mathrm{i} \phi_j \sigma_z} \right] = \left( \begin{array}{cc} P(x) & \mathrm{i} Q(x) \sqrt{1 - x^2}\\ \mathrm{i} Q^*(x) \sqrt{1 - x^2} & P^*(x) \end{array} \right)," /></a>

where <img src="http://chart.googleapis.com/chart?cht=tx&chl= P_{\mathrm{Re}}=f" style="border:none;">.

The package contains two kinds of solvers:

- Optimization-based solver 
- Direct solver (namely the GSLW method and the Haah method)

The package also contains an implementation of the Remez algorithm for finding polynomial approximation.

## Citing our work

If you find our work useful or you use our work in your own project, please consider to cite our work.

- Dong, Y., Meng, X., Whaley, K.B. and Lin, L., 2021. Efficient phase-factor evaluation in quantum signal processing. Physical Review A, 103(4), p.042419.
- Wang, J., Dong, Y. and Lin, L., 2021. On the energy landscape of symmetric quantum signal processing. Quantum 6 (2022): 850.
- Dong, Y., Lin, L., Ni, H., & Wang, J. (2024). Infinite quantum signal processing. Quantum, 8, 1558.
- Dong, Y., Lin, L., Ni, H., & Wang, J. (2024). Robust iterative method for symmetric quantum signal processing in all parameter regimes. SIAM Journal on Scientific Computing, 46(5), A2951-A2971.
- Alexis, M., Lin, L., Mnatsakanyan, G., Thiele, C., & Wang, J. (2024). Infinite quantum signal processing for arbitrary Szeg\H {o} functions. arXiv preprint arXiv:2407.05634.
- Ni, H., Sarkar, R., Ying, L., & Lin, L. (2025). Inverse nonlinear fast Fourier transform on SU (2) with applications to quantum signal processing. arXiv preprint arXiv:2505.12615.

Other references: 

- A. Gilyén, Y. Su, G. H. Low, and N. Wiebe. Quantum singular value transformation and beyond: exponential improvements for quantum matrix arithmetics. In Proceedings of the 51st Annual ACM SIGACT Symposium on Theory of Computing, pages 193–204, 2019
- J. Haah. Product decomposition of periodic functions in quantum signal processing.Quantum, 3:190, 2019



##  Authors

We hope that the package is useful for your application. If you have any bug reports or comments, please feel free to email one of the software authors:

* Yulong Dong, dongyl@berkeley.edu

* Jiasu Wang, jiasu@berkeley.edu

* Xiang Meng, mengxianglgal@gmail.com

* Hongkang Ni, hongkang@stanford.edu

* Lin Lin, linlin@math.berkeley.edu

  

## Installation

- You can download qsppack and run

	`>> startup`

	This adds the folder of the solver to MATLAB's path variable.

- Alternatively you can install qsppack to your current directory by pasting the code below to your MATLAB command window:

    ```matlab
	unzip('https://github.com/qsppack/qsppack/archive/master.zip')
	movefile('QSPPACK-master', 'qsppack')
	addpath(fullfile(cd,'qsppack','Solvers','Optimization')), savepath
    ```

Then you can test qsppack by running

`>> cd Examples`

`>> test_HS`



