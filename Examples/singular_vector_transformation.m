%% Singular vector transformation

%%
% (example/singular_vector_transformation.m)

%% 
% Singular vector transformation implements a matrix that transforms the 
% right singular vector to the corresponding left singular vector. Suppose 
% that we have access to a unitary $U$, its inverse $U^\dagger$ and the 
% controlled reflection operators $(2\Pi-I)$, $(2\tilde{\Pi}-I)$.
% $A:=\tilde{\Pi} U\Pi$ is of our interest and has a singular value decompostion 
% $A=\sum_k \sigma_k |\phi_k \rangle \langle \psi_k|$.
% Singular vector transformation algorithm, first proposed in [GSLW], 
% transforms an arbitrary input state 
% $\sum_k \alpha_k|\psi_k\rangle $ to $ \sum_k \alpha_k|\phi_k\rangle$.

%%
% This algorithm can be used for devising a new method for singular value 
% estimation. It also has many applications, for example, efficient ground 
% state preparation of certain local Hamiltonians.

%%
% As shown in the proof of Theorem 1 in [GSLW], this algorithm aims at
% achieving the singular value transformation of A for sign function $f$, 
% that is, $$f^{(SV)}=\sum_k  |\phi_k\rangle\langle \psi_k|.$$
% In practice, we find an odd polynomial approximation to $f$, denoted as 
% $f_{poly}$, on the interval $D_{\kappa}=[-1, -\delta]\cup [\delta,1]$. 
% By applying the singular value transformation for $f_{poly}$ instead, 
% this algorithm mapps the a right singular vector having singular value at
% least $\delta$ to the corresponding left singular vector.

%%
% For numerical demonstration, we scale the target function by a factor of
% 0.8, to improve the numerical stability.
delta = 0.1;
targ = @(x) 0.8*sign(x);

%%
% We call a subroutine to find the best odd polynomial approximating $f(x)$ 
% on the interval $D_{\kappa}$, where we solves the problem by convex 
% optimization. Here are the parameters set for the subroutine.
opts.intervals=[-1,-delta, delta,1];
opts.objnorm = Inf;
opts.epsil = 0.1;
opts.npts = 500;
opts.isplot= true;
opts.fscale = 1; % disable further rescaling of f(x)

parity = 1;
deg = 251; % agrees with parity  
coef_full=cvx_poly_coef(targ, deg, opts);

%% 
% The solver outputs all Chebyshev coefficients while we have to post-select 
% those of odd order due to the parity constraint.
coef = coef_full(1+parity:2:end);

%% Solving the phase factors
% We use Newton's method for solving phase factors. The parameters of the 
% solver is initiated as follows.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = false;
opts.targetPre = true;
opts.method = 'Newton';
[phi_proc,out] = QSP_solver(coef,parity,opts);

%% Verifying the solution
% We verify the solved phase factors by computing the residual error in 
% terms of the normalized $l^{2}$ norm
% $$residual\_norm = \sqrt{\sum_{k=1,\cdots,K} (g(x_k,\Phi^*)-f_{poly}(x_k))^2}$$
% Using 1000 equally spaced points, the residual error is $ 2.7624e-13$
% which attains almost machine precision. We also plot the pointwise error.

xlist1 = linspace(-1,-delta,500)';
xlist2 = linspace(delta,1,500)';
xlist = cat(1, xlist1,xlist2);
func = @(x) ChebyCoef2Func(x, coef, parity, true);
targ_value = targ(xlist);
func_value = func(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-func_value, 2);
disp('The residual error is');


disp(err);
figure()
plot(xlist,QSP_value-func_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$g(x,\Phi^*)-f_\mathrm{poly}(x)$', 'Interpreter', 'latex')
print(gcf,'singular_vector_transformation_error.png','-dpng','-r500');


