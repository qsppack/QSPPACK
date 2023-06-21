%% Gaussian state preparation

%%
% (example/gaussian_state_preparation.m)

%%
% In this example, We seek to prepare an $N=2^n$ dimensional quantum state
% on $n$ qubits with amplitudes described by function $f(x) =
% e^{-\frac{\beta}{2}x^2}$, which is defined as 
%
% $$ |\Psi\rangle := \frac{1}{\sqrt{\sum_x |f(x)|^2}} \sum_{x=0}^{N-1} f(x) |x\rangle.$$
%
% This is called Gaussian states and are widely used in quantum algorithms,
% e.g. in chemistry, simulation of quantum field theories, and finance.

%%
% For numerical demonstration, we consider target function to be 
% $f(x) = e^{-50 x^2}$.
%
% According to the paper ``Quantum state preparation without coherent
% arithmetic", we first apply the circuit $U_{sin}$ that block-encodes
%
% $$ \sum_x \sin(x/N)|x\rangle\langle x|.$$
%
% In order to achieve theb desired 
% transformation $f$, we need to find the phase factors for 
% $h(z) = f(\sin^{-1}(z))$ instead. Here, we introduce a new  variable 
% $z = \sin(x)$ and the domain is $z\in [0,\sin(1)]$.

beta = 100;
targ = @(x) exp(-beta/2 *asin(x).^2);


%%
% To numerically find the best even polynomial approximating $h(z)$,  we use
% a subroutine which solves the problem using convex optimization. Here are
% the parameters of the subroutine.

deg = 100;
opts.intervals=[0,sin(1)];
opts.objnorm = Inf;
opts.epsil = 0.01;
opts.npts = 500;
opts.fscale = 0.99;
opts.isplot=true;

%%
% This subroutine yields the coefficients of the approximation polynomial 
% in the Chebyshev basis. we have to post-select those of odd order due 
% to the parity constraint.
coef_full=cvx_poly_coef(targ, deg, opts);
parity = mod(deg, 2);
coef = coef_full(1+parity:2:end);


%% Solving the phase factors
% We use LBFGS method for solving phase factors. The parameters of the 
% solver is initiated as follows.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = true;
opts.targetPre = true;
opts.method = 'LBFGS';

[phi_proc,out] = QSP_solver(coef,parity,opts);

%% Verifying the solution
% We verify the solved phase factors by computing the residual error in 
% terms of $l^{\infty}$ norm.

xlist = linspace(0,sin(1),500)';
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
print(gcf,'gaussian_state_preparation.png','-dpng','-r500');


