%% Quantum Gibbs state preparation

%%
% (example/quantum_gibbs_state_preparation.m)

%% 
% The preparation of a quantum Gibbs state is an essential part of quantum 
% computation and has wide-ranging applications in various areas, including 
% quantum simulation, quantum optimization, and quantum machine learning. 
% The Gibbs state for a quantum Hamiltonian $H$ is defined as the density 
% operator 
%
% $$\rho = \frac{e^{-\beta H}}{tr(e^{-\beta H})}$$
%
% The basic idea for preparing Gibbs state is approximately implementing 
% the map $H \to e^{-\beta H}$. 

%%
% In QSP, we only need to consider the scalar function $f(x)=e^{-\beta x}$. 
% Without loss of generality, we assume the matrix is normalized so that 
% $\|H\|\leq 1$. Besides, the eigenvalues of $H$ lie in the interval
% $D_{\delta}:=[\delta, 1]$. 
% Hence, we need to find an odd polynomial approximation over $D_{\delta}$.
% In our numerical example, we choose $\beta =2$ and $\delta=0.1$.

beta = 2;
targ = @(x) exp(-beta * x);
delta = 0.2;
deg = 151;

%%
% To numerically find the best polynomial approximating $f(x)$ in the sense
% of 2 norm on the interval $D_{\delta}$, we use a subroutine which solves
% the problem using convex optimization. We first set the parameters of the
% subroutine.

opts.intervals=[delta,1];
opts.objnorm = 2;
% opts.epsil is usually chosen such that target function is bounded by 1-opts.epsil over D_delta 
opts.epsil = 0.2;
opts.npts = 500;
opts.fscale = 1;
opts.isplot = true;
%%
% This subroutine yields the coefficients of the approximation polynomial 
% in the Chebyshev basis.
coef_full=cvx_poly_coef(targ, deg, opts);

%%
% The subroutine outputs all Chebyshev coefficients while we have to 
% post-select those of odd order due to the parity constraint.
parity = mod(deg, 2);
coef = coef_full(1+parity:2:end);

%%
% To visualize this polynomial approximation, we plot it with the target function.

% convert the Chebyshev coefficients to Chebyshev polynomial
func = @(x) ChebyCoef2Func(x, coef, parity, true);
xlist1 = linspace(0.5*delta,1,500)';
func_value1 = func(xlist1);
targ_value1 = targ(xlist1);
xlist2 = linspace(-1,1,1000)';
func_value2 = func(xlist2);


%% Solving the phase factors
% We use Newton's method for solving phase factors. The parameters of the 
% solver is initiated as follows.

opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = false;
opts.targetPre = true;
opts.method = 'Newton';
% solve phase factors
[phi_proc,out] = QSP_solver(coef,parity,opts);


%% Verifying the solution
% We verify the solved phase factors by computing the residual error in 
% terms of $l^{\infty}$ norm.

xlist = linspace(delta,1,500)';
func = @(x) ChebyCoef2Func(x, coef, parity, true);
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


