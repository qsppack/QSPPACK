%% Matrix inversion

%%
% (example/matrix_inversion.m)

%% 
% In matrix inversion, the function of interest is $f(x)=1/x$. To
% implement QSP, we need a polynomial approximation of $1/x$ over interval
% $D_{\kappa}:=[1/\kappa, 1]$. Here $\kappa>1$ is the condition number
% of a matrix.

%%
% As an example, consider $\frac{1}{20x}$
kappa = 10;
targ = @(x) (1/(2*kappa))./x;

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
opts.epsil = 0.1;
opts.npts = 500;
% opts.isplot = false;
deg = 121;

%%
% Since the $L^{\infty}$ norm of target function is bounded by 0.5, 
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

figure()
plot(xlist,QSP_value-func_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$g(x,\Phi^*)-f_\mathrm{poly}(x)$', 'Interpreter', 'latex')
print(gcf,'quantum_linear_system_problem_error.png','-dpng','-r500');

%%
% Show the quality of polynomial approximation.
figure()
hold on
xlist1 = linspace(1/kappa,1,500)';
targ_value1 = targ(xlist1);
plot(xlist1,opts.fscale * targ_value1,'b-','linewidth',2)
plot(-xlist1,-opts.fscale * targ_value1,'b-','linewidth',2)
xlist1 = linspace(-1,1,1000)';
func_value1 = func(xlist1);
plot(xlist1,func_value1,'-.')
hold off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$f(x)$', 'Interpreter', 'latex')
legend('target', '', 'polynomial',...
  'location','se')

figure()
plot(xlist,func_value-opts.fscale * targ_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$f_\mathrm{poly}(x)-f(x)$', 'Interpreter', 'latex')
print(gcf,'quantum_linear_system_problem_polynomial.png','-dpng','-r500');




