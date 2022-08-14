%% Quantum Gaussian filter
% Jiasu Wang, July 2022

%%
% (example/quantum_gaussian_filter.m)

%% 
% The Quantum Gaussian filter solves the approximate ground state by 
% performing a Gaussian function of Hamiltonian $H$ on a given initial 
% state which has a sufficient overlapping with the ground state. 
% The Gaussian filter operator is $e^{-(H-\mu I)^2/\sigma^2}$, where the
% expected value $\mu$ and variance $\sigma^2$ of the Gaussian function 
% correspond to minus shift-energy and width of the Gaussian filter.
% Let $|\lambda_j\rangle$ be the eigenstate of $H$ corresponding to j-th 
% smallest eigenvalue $\lambda_j\rangle$. The Gaussian filter operator 
% results in an additional weight $e^{-(\lambda_i-\mu I)^2/\sigma^2}$ for 
% each eigenstate, which monotonically decreases with the eigenvalue. 
% Without loss of generlaity, we assume that all eigenvalues lie in $[0,1]$.
% Then as the variance $\sigma^2$ or the expected value $\mu$ of the 
% Gaussian function decreases, the resulting state converges to the ground 
% state $|\lambda_0\rangle$.

%%  
% For numerical demonstration, we consider even function 
%
% $$0.8*e^{-(|x|-0.5)^2/3}$$ 
%
% whose $l^{\infty}$ norm over $[-1,1]$ is strictly bounded by $0.8$.

parity = 0;
targ = @(x) 0.8*exp(-(abs(x)-0.5).^2/3);
% Compute its Chebyshev coefficients. 
d = 200;
f = chebfun(targ,d);
coef = chebcoeffs(f);
% Only need part of Chebyshev coefficients.
coef = coef(parity+1:2:end);

%% Solving the phase factors
% We use LBFGS method for solving phase factors. The parameters of the 
% solver is initiated as follows.

opts.maxiter = 100;
opts.criteria = 1e-12;
opts.targetPre = true;
opts.method = 'Newton';
% use the real representation to speed up the computation
opts.useReal = true;
% solve phase factors
[phi_proc,out] = QSP_solver(coef,parity,opts);

%% Verifying the solution
% We verify the solved phase factors by computing the residual error in 
% terms of the normalized $l^{\infty}$ norm
% $$residual\_norm = \max_{k=1,\cdots,K} |g(x_k,\Phi^*)-f_{poly}(x_k)|$$
% Using 1000 equally spaced points, the residual error is $5.2558e-13$
% which attains almost machine precision. We also plot the pointwise error.

xlist = linspace(0,1,1000)';
func = @(x) ChebyCoef2Func(x, coef, parity, true);
targ_value = targ(xlist);
func_value = func(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-func_value,Inf);
disp('The residual error is');
disp(err);


figure()
plot(xlist,QSP_value-func_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$g(x,\Phi^*)-f_\mathrm{poly}(x)$', 'Interpreter', 'latex')
print(gcf,'quantum_gaussian_filter_error.png','-dpng','-r500');



