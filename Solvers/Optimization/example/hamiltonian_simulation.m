%% Hamiltonian simulation
% Jiasu Wang, July 2022

%%
% (example/hamiltonian_simulation.m)

%% Hamiltonian simulation
% In Hamiltonian simulation, the function of our interest is $f(x)=e^{-i\tau x}$.
% In practice, one always consider implementing the real and imaginary component of the complex polynomial individually
% and then combine them by linear combination of unitaries.
% Thus we only need to determine the phase factors corresponding to those polynomials approximating $\cos(\tau x)$ and $\sin(\tau x)$. 

%% Approxiomating the real compenent
% As an example, consider function $0.5*\cos(100 x)$, whose $L^{\infty}
% norm over $[-1,1]$ is strictly bounded by $\frac{1}{2}$.
parity = 0;
tau = 100;
targ = @(x) 0.5*cos(tau.*x);

%%
% |chebfun| provides its Chebyshev coefficients. We truncate the series up
% to $d=1.4| \tau |+\log(1/\epsilon_0)$ such that the approximation error is 
% bounded by $\epsilon_0$.
d = ceil(1.4*tau+log(1e14));
f = chebfun(targ,d);
coef = chebcoeffs(f);

%%
% We only need its Chebyshev coefficients with respect to $T_{2k}$, where
% $k$ is nonegative integer.
coef = coef(parity+1:2:end);

%%
% Set up the parameters for solver.
opts.maxiter = 100;
opts.criteria = 1e-12;

%%
% Set |opts.useReal| to be |true| will increase the computing speed.
opts.useReal = false;

%%
% We want the real component of the upper left entry of the QSP unitary
% matrix to be the target function.
opts.targetPre = true;

%%
% Use Contraction mapping method to find phase factors
opts.method = 'CM';
[phi_proc,out] = QSP_solver(coef,parity,opts);

%%
% We do the following test to demonstrate that the obtained phase factors 
% satisfy expectation.
xlist = rand(100,1)*2-1;
targ_value = targ(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-targ_value,1)/100;
disp('the residual error is');
disp(err);

