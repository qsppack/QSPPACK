%% Negative power function
% Jiasu Wang, August 2022

%%
% (example/negative_power_function.m)

%% 
% As a generalization of the quantum linear system problem (QLSP), we 
% consider the following problem: given the access to an invertible matrix
% $A$ and a normalized quantum state $|b\rangle$, we want to construct a
% quantum state 
%
% $$|x\rangle:=\frac{A^{-c}|b\rangle}{\|A^{-c}|b\rangle\|_2},$$
%
% where $c$ is a interger larger than 1. Similarly, the implement boils
% down to a scalar function $f(x)=x^{-c}$. Without loss of generality, we 
% assume the matrix is normalized. We have to find a polynomial with parity
% approximating $f(x)$ over interval
%
% $$D_{\kappa}:=[-1,-\kappa]\cup[-\kappa,1].$$

%%
% For numerical demonstration, we set $\kappa =10$, $c=2$ and scale down the 
% target function by a factor of $\frac{1}{2\times 10^4}$ to meet the
% $l^{\infty}$ norm condition over $[-1,1]$.
%
% [LL: why the normalization factor is so large here? Should it just be]
% bounded in $D_{\kappa}$?]
kappa = 10;
targ = @(x) (1/(2*kappa^2))/x.^2;
deg = 150; % approximate f(x) by a polynomial of degree deg
parity = 0;

%%
% We want to find the best approximation polynomial with degree up to $d$ 
% in terms of $L^{1}$ norm. To achieve it, we call |cvx_poly_coef|, 
% which solves the optimization problem and outputs the Chebyshev
% coefficients of the best approximation polynomial.

%%
%
% $$\min_{f\in R[x], \deg(f)\leq d} \max_{x\in [0,1]} |f(x)-sign(x-a)|$$
%
% subject to $\max_{x\in[0,1]} |f(x)|\leq 1-epsil$.
opts.intervals=[1/kappa,1];
opts.objnorm = Inf;
opts.epsil = 0.2;
opts.npts = 400;
opts.isplot = false;
opts.fscale = 1;

%%
% The inf norm of target function exceeds 1. Hence we need rescale it.
opts.fscale = 5e-5;

%%
% Call |cvx_poly_coef| to compute the Chebyshev coefficients for the best
% approximation even polynomials
coef_full=cvx_poly_coef(targ, deg, opts);

%% 
% We only need part of its Chebyshev coefficients.
coef = coef_full(1+parity:2:end);

%%
% Set up the parameters for solver.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = false;
opts.targetPre = true;

%%
% Use Newton method to find phase factors
opts.method = 'Newton';
[phi_proc,out] = QSP_solver(coef,parity,opts);

%% Verifying the solution
%
% We verify the solved phase factors by computing the residual error in 
% terms of the normalized $l^{2}$ norm
%
% $$\mathrm{res} =\sqrt{ \sum_{k=1,\cdots,K} (g(x_k,\Phi^*)-f_{poly}(x_k))^2}$$
%
% Using 1000 equally spaced points, the residual error is $5.2558e-13$
% which attains almost machine precision. We also plot the pointwise error.
xlist1 = linspace(-1,-1/kappa,500)';
xlist2 = linspace(1/kappa,1,500)';
xlist = cat(1, xlist1,xlist2);
func = @(x) ChebyCoef2Func(x, coef, parity, true);
targ_value = targ(xlist);
func_value = func(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-func_value,2);
disp('The residual error is');
disp(err);


figure()
plot(xlist,QSP_value-func_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$g(x,\Phi^*)-f_\mathrm{poly}(x)$', 'Interpreter', 'latex')
print(gcf,'negative_power_function_error.png','-dpng','-r500');

