%% Generate the phase factors for a Chebyshev polynomial
% Phase factor for a simple 3rd order Chebyshev polynomial
%
% $f(x) = 0.2 T_1(x) + 0.4 T_3(x)$

%%
% (example/chebyshev_ex2.m)

deg = 3;
parity = mod(deg,2);
coef_targ = zeros(deg+1,1);
coef_targ(2)=0.2;
coef_targ(4)=0.4;
targ = chebfun(coef_targ, 'coeffs');

%%
% We only need its Chebyshev coefficients with respect to $T_{2k}$, where
% $k$ is nonegative integer.
coef = coef_targ(parity+1:2:end);

%%
% Set up the parameters for the solver.
opts.maxiter = 100;
opts.criteria = 1e-12;

%%
% Set |opts.useReal| to be |true| will increase the computing speed.
opts.useReal = true;

%%
% We want the real component of the upper left entry of the QSP unitary
% matrix to be the target function.
opts.targetPre = true;

%%
% Use the fixed point iteration method to find phase factors
opts.method = 'Newton';
[phi_proc,out] = QSP_solver(coef,parity,opts);

disp('Symmetric phase factors = ');
disp(phi_proc);


%%
% We do the following test to demonstrate that the obtained phase factors 
% satisfy expectation.
xlist = linspace(-1, 1, 1000)';
targ_value = targ(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-targ_value,1)/length(xlist);
disp('The residual error is');
disp(err);

figure
hold on
plot(xlist,QSP_value,'b-')
plot(xlist,targ_value,'r--')
hold off
legend('QSP','Target')
xlabel('$x$', 'Interpreter', 'latex')

figure
plot(xlist,QSP_value-targ_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$g(x,\Phi^*)-f(x)$', 'Interpreter', 'latex')


