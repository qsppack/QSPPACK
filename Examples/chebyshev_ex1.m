%% Generate the phase factors for a Chebyshev polynomial
% Phase factor for Chebyshev polynomials of the first kind is a zero
% vector.

%%
% (example/chebyshev_ex1.m)

k = 5;
deg = 2*k;
parity = mod(deg,2);
coef_targ = zeros(deg+1,1);
coef_targ(deg+1)=1;
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
xlist = linspace(0, 1, 1000)';
targ_value = targ(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-targ_value,1)/length(xlist);
disp('The residual error is');
disp(err);

plot(xlist,QSP_value-targ_value)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$g(x,\Phi^*)-f(x)$', 'Interpreter', 'latex')
print(gcf,'hamiltonian_simulation.png','-dpng','-r500');

