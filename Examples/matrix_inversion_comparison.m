%% Matrix Inversion Comparison
% This script compares the performance of direct construction and convex optimization-based methods
% for approximating 1/x (properly scaled). The direct construction
% follows the recipe from

%% Initialization
kappa = 5; % Condition number
beta = 2;
x = chebfun('x');
targ = (1/(beta*kappa))/x;

%%
% We want to find the best approximation polynomial with degree up to $d$ 
% in terms of $L^{\infty}$ norm. To achieve it, we call |cvx_poly_coef|, 
% which solves the optimization problem and outputs the Chebyshev
% coefficients of the best approximation polynomial.

%%
%
% $$\min_{f\in R[x], \deg(f)\leq d} \max_{x\in D_{\kappa}} |f(x)-1/x|$$
%
% subject to $\max_{x\in[0,1]} |f(x)|\leq 1-\epsilon$.
opts.intervals=[1/kappa,1];
opts.objnorm = Inf;
opts.epsil = 0.01;
opts.npts = 500;
opts.isplot = true;
deg = 51;

%%
% Since the $L^{\infty}$ norm of target function is bounded by 0.5, 
% we don't need rescale the target function and |opts.fscale| is set to be
% 1.
opts.fscale = 1;
coef_full=cvx_poly_coef(targ, deg, opts);
p_cvx = chebfun(coef_full, 'coeffs');

%% 
% We only need part of its Chebyshev coefficients.
parity = mod(deg, 2);
coef = coef_full(1+parity:2:end);

%%
% Set up the parameters for solver.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = true;
opts.targetPre = true;

%%
% Use Newton method to find phase factors
opts.method = 'Newton';
[phi_proc,out] = QSP_solver(coef,parity,opts);

%%
% We do the following test to demonstrate that the obtained phase factors 
% satisfy expectation.
xlist = linspace(1/kappa,1,1000)';
func = @(x) ChebyCoef2Func(x, coef, parity, true);
targ_value = targ(xlist);
func_value = func(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-func_value,Inf);
disp('The residual error is');
disp(err);


%% Compare with the direct construction algorithm

deg_direct = 81;
% epsilon_direct = 1e-6; % Desired accuracy
% b = ceil(kappa^2 * log(kappa/epsilon_direct));
b = 200;

if(1)
  f_approx = (1-(1-x^2)^b)/x;
  coef_direct = chebcoeffs(1/(beta*kappa)*f_approx, deg_direct);
else
  % Remark: this Chebyshev expansion formula is not numerically stable and
  % should not be used
  d = ceil(sqrt(b * log(4*b/epsilon_direct)));
  coef_direct = zeros(deg_direct+1,1);
  for j = 0 : (deg_direct-1)/2
    sum_coef = 0;
    for i = j+1:b
      sum_coef = sum_coef + nchoosek(2*b, b+i);
    end
    coef_direct(2*(j+1)) = 4/(beta*kappa) * ((-1)^j * sum_coef / 2^(2*b));
  end
end

p_direct = chebfun(coef_direct, 'coeffs');

str_direct = sprintf('Analytic   Deg = %d', deg_direct);
str_cvx    = sprintf('Optimized  Deg = %d', deg);


figure;
plot(x, p_direct, 'r-.', x, p_cvx, 'b-', x, targ, 'k--');
title(sprintf('$\\beta$=%d, $\\kappa$=%d', beta, kappa));
legend(str_direct, str_cvx, 'Target $(\beta\kappa x)^{-1}$');
xlabel('$x$');
xlim([0.01, 1]);
ylim([0,1]);
print(sprintf('matrix_inversion_comparison.png', deg),'-dpng','-r500');

figure
plot(x, abs(p_direct - targ), 'r-.', x, abs(p_cvx - targ), 'b-');
legend(str_direct, str_cvx);
title(sprintf('Error $\\kappa$=%d', kappa));
xlabel('$x$');
xlim([1/kappa, 1]);
ylim([0, 1e-4]);
print(sprintf('matrix_inversion_comparison_error.png', deg),'-dpng','-r500');


%% References
% [1] Childs, A. M., Kothari, R., & Somma, R. D. (2017). Quantum
% Algorithm for Systems of Linear Equations with Exponentially Improved
% Dependence on Precision. SIAM Journal on Computing, 46(6), 1920–1950.
% https://doi.org/10.1137/16M1087072
