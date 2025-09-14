%% Hamiltonian simulation

%%
% (example/hamiltonian_simulation.m)

%% Hamiltonian simulation
% In Hamiltonian simulation, the function of interest is $f(x)=e^{-i\tau
% x}$.  In practice, the real and imaginary component of the complex
% polynomial is implemented separately, and are then combined by linear
% combination of unitaries.
% Thus we only need to determine the phase factors corresponding to those polynomials approximating $\cos(\tau x)$ and $\sin(\tau x)$. 

%% Approxiomating the real compenent
% Consider the real part $0.5\cos(100 x)$, whose $L^{\infty}$
% norm over $[-1,1]$ is strictly bounded by $\frac{1}{2}$.
parity = 0;
tau = 100;
targ = @(x) 0.5*cos(tau.*x);

%%
% The Chebyshev coefficients can be computed using |chebfun|. We truncate the series up
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
opts.method = 'FPI';
%opts.method = 'Newton';

[phi_proc,out] = QSP_solver(coef,parity,opts);

%%
% We do the following test to demonstrate that the obtained phase factors 
% satisfy expectation.
xlist = linspace(0, 1, 1000)';
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
print(gcf,'hamiltonian_simulation.png','-dpng','-r500');

%%
% Plot the phase factor
figure
plot(phi_proc,'b-o')
xlabel('$i$', 'Interpreter', 'latex')
ylabel('$\phi_i$', 'Interpreter', 'latex')
title('Phase factor');

%%
% Demonstrate decay behavior
figure
phi_shift = phi_proc(1:end);
phi_shift(1) = phi_shift(1)-pi/4;
phi_shift(end) = phi_shift(end)-pi/4;
semilogy(abs(phi_shift))
axis tight
ylabel('$|\Phi-\Phi_0|$')
title('Decay behavior')

%%
% Generate the required figure for the LaTeX document with 3 subfigures
figure('Position', [100, 100, 1200, 400])

% Left subplot: polynomial approximation vs actual function
subplot(1,3,1)
hold on
plot(xlist, targ_value, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Target')
plot(xlist, QSP_value, 'b--', 'LineWidth', 1.5, 'DisplayName', 'QSP')
hold off
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 12)
ylim([-1,1])
legend('Location', 'best')
grid on
box on

% Middle subplot: error between polynomial and actual function
subplot(1,3,2)
plot(xlist, QSP_value - targ_value, 'k-', 'LineWidth', 1.5)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Error', 'FontSize', 12)
grid on
box on

% Right subplot: phase factors (after removing pi/4 factor)
subplot(1,3,3)
plot(1:length(phi_shift), phi_shift, 'bo-', 'MarkerSize', 4, 'LineWidth', 1)
xlabel('Index $j$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Phase Factors', 'FontSize', 12)
grid on
axis tight
box on

% Save the figure
print(gcf, 'qsp_cos100x.png', '-dpng', '-r300')
