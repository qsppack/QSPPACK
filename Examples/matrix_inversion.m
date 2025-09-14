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
opts.isplot = true;
deg = 101;

display(['Degree of the approximating polynomial: ', num2str(deg)]);

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
% Generate the required figure for the LaTeX document with 3 subfigures
% Create extended x-list from 0 to 1 for QSP plotting
xlist_extend = linspace(0, 1, 1000)';
QSP_value_extend = QSPGetEntry(xlist_extend, phi_proc, out);

figure('Position', [100, 100, 1200, 400])

% Left subplot: target function on original interval and QSP on extended interval
subplot(1,3,1)
hold on
plot(xlist, targ_value, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Target')
plot(xlist_extend, QSP_value_extend, 'b--', 'LineWidth', 1.5, 'DisplayName', 'QSP')
hold off
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 12)
ylim([0,1])
legend('Location', 'best')
grid on

% Middle subplot: difference between QSP and target function
subplot(1,3,2)
plot(xlist, QSP_value - func_value, 'k-', 'LineWidth', 1.5)
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('QSP error', 'FontSize', 12)
grid on
% subplot(1,3,2)
% plot(xlist, QSP_value - targ_value, 'k-', 'LineWidth', 1.5)
% xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12)
% ylabel('QSP - Target', 'FontSize', 12)
% title('QSP vs Target Error', 'FontSize', 12)
% grid on

% Right subplot: phase factors (after removing pi/4 factor)
subplot(1,3,3)
phi_shift = phi_proc(1:end);
phi_shift(1) = phi_shift(1)-pi/4;
phi_shift(end) = phi_shift(end)-pi/4;
semilogy(1:length(phi_shift), abs(phi_shift), 'bo-', 'MarkerSize', 4, 'LineWidth', 1)
xlabel('Index $j$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$|\psi_j|$', 'Interpreter', 'latex', 'FontSize', 12)
grid on
axis tight
box on

% Save the figure
print(gcf, 'qsp_inversion.png', '-dpng', '-r300')

