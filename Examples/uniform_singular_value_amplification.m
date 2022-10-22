%% Uniform singular value amplification

%%
% (example/uniform_singular_value_amplification.m)

%% 
% Uniform singular value amplification uniformly amplify the singular 
% values of a matrix represented as a projected unitary.
% Suppose that we have access to a unitary $U$, its inverse $U^\dagger$ and the 
% controlled reflection operators $(2\Pi-I)$, $(2\tilde{\Pi}-I)$. 
% $A:=\tilde{\Pi} U\Pi$ is of our interest and has a singular value decompostion 
% $A=W\Sigma V$. The goal of uniform singular value amplification is to
% implement the singular value transformation of $A$ for function $f(x)$.
% Here, $f(x)$ is linear function $\lambda x$ over some interval and takes 
% value zero for other x. 

%%
% For numerical demonstration, we consider target function $f(x)=\lambda
% \max(x-a,0)$, where $\lambda =1.5$ and $a=0.5$. 
a = 0.5;
lambda = 1.5;
targ = @(x) lambda*max(abs(x)-a,0);

%%
% To numerically find the best even polynomial approximating $f(x)$ in the sense
% of 2 norm on the interval $D_{\delta}$, we use a subroutine which solves
% the problem using convex optimization. We first set the parameters of the
% subroutine.

deg = 500;
delta =0.01;
opts.intervals=[0,a-delta, a+delta,1];
opts.objnorm = Inf;
% opts.epsil is usually chosen such that target function is bounded by 1-opts.epsil over D_delta 
opts.epsil = 0.2;
opts.npts = 500;
opts.fscale = 1;

%%
% This subroutine yields the coefficients of the approximation polynomial 
% in the Chebyshev basis.
coef_full=cvx_poly_coef(targ, deg, opts);

%%
% The subroutine outputs all Chebyshev coefficients while we have to 
% post-select those of odd order due to the parity constraint.
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
% terms of the normalized $l^{\infty}$ norm
% $$residual\_norm = \max_{k=1,\cdots,K} |g(x_k,\Phi^*)-f_{poly}(x_k)|$$
% Using 1000 equally spaced points, the residual error is $1.3767e-14$
% which attains almost machine precision. We also plot the pointwise error.
xlist1 = linspace(0,a-delta,500)';
xlist2 = linspace(a+delta,1,500)';
xlist = cat(1, xlist1,xlist2);
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


