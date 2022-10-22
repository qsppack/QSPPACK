%% Quantum Gaussian filter
% Jiasu Wang, July 2022

%%
% (example/quantum_gaussian_filter.m)

%% 
% The Quantum Gaussian filter can be used to approximate the ground state by 
% performing a Gaussian function to a Hamiltonian $H$ on a given initial 
% state. The Gaussian filter operator is $e^{-(H-\mu I)^2/\sigma^2}$.

%%  
% For numerical demonstration, we consider even function 
%
% $$0.99 * e^{-(|x|-0.5)^2/0.1^2}$$ 
%
% whose $l^{\infty}$ norm over $[-1,1]$ is strictly bounded by $0.99$, and
% this is very close to the fully coherent regime.

parity = 0;
targ = @(x) 0.99 * exp(-(abs(x)-0.5).^2/0.1^2);
% Compute its Chebyshev coefficients. 
d = 100;
f = chebfun(targ,d);
coef = chebcoeffs(f);
% Only need part of Chebyshev coefficients.
coef = coef(parity+1:2:end);

%% Solving the phase factors
% We use the Newton method for solving phase factors. The parameters of the 
% solver is initiated as follows.

opts.maxiter = 100;
opts.criteria = 1e-12;
opts.targetPre = true;
opts.method = 'Newton';
% use the real representation to speed up the computation
opts.useReal = true;
% solve phase factors
[phi_proc,out] = QSP_solver(coef,parity,opts);

%% Verifying the solution
% We verify the solved phase factors by computing the residual error in 
% terms of the normalized $l^{\infty}$ norm.

xlist = linspace(0,1,1000)';
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
print(gcf,'quantum_gaussian_filter_error.png','-dpng','-r500');

%%
% Show the quality of polynomial approximation.
figure()
hold on
targ_value = targ(xlist);
plot(xlist,targ_value,'b-','linewidth',2)
plot(xlist,func_value,'-.')
hold off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$f(x)$', 'Interpreter', 'latex')
legend('target',  'polynomial',...
  'location','se')

figure()
plot(xlist,func_value-targ_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$f_\mathrm{poly}(x)-f(x)$', 'Interpreter', 'latex')
print(gcf,'quantum_linear_system_problem_polynomial.png','-dpng','-r500');


