%% Matrix inversion
% Jiasu Wang, July 2022

%%
% (example/matrix_inversion.m)

%% 
% In Matrix inversion, the function of our interest is $f(x)=1/x$. To
% implement QSP, we need a polynomial approximation of $1/x$ over interval
% $D_{\kappa}:=[1/\kappa, 1]$.

%%
% As an example, consider function $\frac{1}{20x}$
targ = @(x) 0.05*1./x;

%%
% We want to find the best approximation polynomial with degree up to $d$ 
% in terms of $L^{\infty}$ norm. To achieve it, we call |cvx_poly_coef|, 
% which solves the optimization problem and outputs the Chebyshev
% coefficients of the best approximation polynomial.

%%
% $\min_{f\in R[x], deg(f)\leq d} \max_{x\in D_{\kappa}} |f(x)-1/x|$
% subject to $\max_{x\in[0,1]} |f(x)|\leq 1-\epsilon$.
kappa = 10;
opts.intervals=[1/kappa,1];
opts.objnorm = Inf;
opts.epsil = 0.2;
opts.npts = 500;
opts.isplot = false;
deg = 150;

%%
% Since the $L^{\infty}$ norm of target function is bounded by 0.5, hence
% we don't need rescale the target function and |opts.fscale| is set to be
% 1.
opts.fscale = 1;
coef_full=cvx_poly_coef(targ, deg, opts);

%% 
% We only need part of its Chebyshev coefficients.
parity = mod(deg, 2);
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

%%
% We do the following test to demonstrate that the obtained phase factors 
% satisfy expectation.
xlist = rand(10,1)*0.9+0.1;
func = @(x) ChebyCoef2Func(x, coef, parity, true);
func_value = func(xlist);
targ_value = targ(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-targ_value,Inf);
disp('the residual error is');
disp(err);

