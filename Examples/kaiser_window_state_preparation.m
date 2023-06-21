%% Kaiser window state preparation

%%
% (example/kaiser_window_state_preparation.m)

%%
% In this example, We seek to prepare the Kaiser window state. Ihe state is
% defined as 
%
% $$ |\Psi\rangle := \frac{1}{\sqrt{\sum_x |f(x)|^2}} \sum_{x=0}^{N-1} f(x) |x\rangle.$$
%
% and the amplitude is described by the Kaiser window function 
% $f(x) = \frac{I_0(\beta \sqrt{1-x^2})}{I_0 (\beta)}$. Here, $I_0$ is the
% zeroth modified Bessel function of the first kind.
% The Kaiser window function can be used in quantum phase estimation to
% boost the success probability without having to (coherently) calculate
% the median of several phase evaluations.


%%
% For numerical demonstration, we consider target function to be 
%
% $$f(x) = \frac{I_0 (8\sqrt(1-2x^2))}{I_0 (8)}.$$
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

beta = 8;
targ = @(x) besselj(0,1j*beta*sqrt(1-2*asin(x).^2))/besselj(0,1j*beta);


%%
% To numerically find the best even polynomial approximating $h(z)$,  we use
% a subroutine which solves the problem using convex optimization. We first
% set the parameters of the subroutine.

deg = 100;
delta =0.01;
opts.intervals=[0,sin(1)];
opts.objnorm = Inf;
opts.epsil = 0.01;
opts.npts = 500;
opts.fscale = 0.98; % scaling factor for the infinity norm
opts.isplot=true;

%%
% This subroutine yields the coefficients of the approximation polynomial 
% in the Chebyshev basis. We post-select those of odd order due to the
% parity constraint

coef_full=cvx_poly_coef(targ, deg, opts);
parity = mod(deg, 2);
coef = coef_full(1+parity:2:end);

%% Solving the phase factors
% We use Newton method for solving phase factors. The parameters of the 
% solver is initiated as follows.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = true;
opts.targetPre = true;
opts.method = 'Newton';

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
print(gcf,'uniform_singular_value_amplification_error.png','-dpng','-r500');


