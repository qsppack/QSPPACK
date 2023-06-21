%% Singular value threshold projectors

%%
% (example/singular_value_threshold_projectors.m)

%%
% Siungular value threshold projectors project out singular vectors with 
% singular values below/above a certain threshold. Suppose 
% that we have access to a unitary $U$, its inverse $U^\dagger$ and the 
% controlled reflection operators $(2\Pi-I)$, $(2\tilde{\Pi}-I)$. 
% $A:=\tilde{\Pi} U\Pi$ is of our interest and has a singular value decompostion 
% $A=W\Sigma V$. For $S\subset R$ let $\Sigma_S$ be the matrix obtained
% from $\Sigma$ by replacing all diagonal entries $\Sigma_{ii}\in S$ by 1
% and all diagonal entries $\Sigma_{ii}\notin S$ by 0. We define
% $\Pi_S:=\Pi V\Sigma_S V^{\dagger}\Pi$ and similarly 
% $\tilde{\Pi}_S:=\tilde{\Pi} W\Sigma_S W^{\dagger}\tilde{\Pi}$.
%
% For example, set $S=[0,\delta]$ and $\Pi_S$ projects out right
% singular vectors with singular value at most $\delta$.
% These threshold projectors play a major role in quantum algorithms.

%% 
% These threshold projectors play a major role in quantum algorithms. One 
% of the applications is the singular value discrimination problem, which 
% decides whether a given quantum state has singular value at most $a$ or 
% it at least $b$.

%%
% As an example, we implement $\Pi_{[0,0.5]}$. We need 
% to implement singular value transformation of $A$ for the following 
% rectangle function. 
targ = @(x) rectangularPulse(-0.5,0.5,x);

%%
% We call a subroutine to find the best even polynomial approximating
% $f(x)$ on the interval
%
% $$D_{\delta}=[0, 0.5-\delta]\cup [0.5+\delta,1].$$
%
% We solve the problem by convex optimization. Here are the parameters set
% for the subroutine.
delta=0.05;
opts.intervals=[0,0.5-delta,0.5+delta,1];
opts.objnorm = Inf;
opts.epsil = 0.1;
opts.npts = 500;
opts.isplot= true;
opts.fscale = 0.9; % scaling factor for the infinity norm

parity = 0;
deg = 250;
coef_full=cvx_poly_coef(targ, deg, opts);

% To meet the $l^{\infty}$ norm requirement, the solver may scale the target
% function and the scaling factor is saved in parameter opts.scale. 
% scaled_targ = @(x) opts.fscale * exp(-beta * x);

%% 
% The solver outputs all Chebyshev coefficients while we have to post-select 
% those of odd order due to the parity constraint.
coef = coef_full(1+parity:2:end);

%% Solving the phase factors
% We use LBFGS method for solving phase factors. The parameters of the 
% solver is initiated as follows.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = false;
opts.targetPre = true;
opts.method = 'Newton';

[phi_proc,out] = QSP_solver(coef,parity,opts);

%% Verifying the solution
% We verify the solved phase factors by computing the residual error in 
% terms of $l^{\infty}$ norm.



xlist1 = linspace(0,0.5-delta,500)';
xlist2 = linspace(0.5+delta,1,500)';
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
print(gcf,'singular_value_threshold_projectors_error.png','-dpng','-r500');
