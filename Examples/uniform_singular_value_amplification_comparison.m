%% Uniform Singular Value Amplification Comparison
% This script compares the performance of direct construction and convex optimization-based methods
% for uniform singular value amplification. It aims to approximate a linear function within a specified
% interval using these methods and evaluates their performance.

%% Initialization
fscale = 1./(1+1.05);       % Scaling factor to ensure bounds
a = 0.2;            % Upper limit of the interval [0, a]
parity = 1;         % Ensures the function is odd
targ = @(x) x/a;    % Target linear function scaled by 1/a

%% Analytic Method using Odd Polynomial Approximation
% Approximates the truncated linear function on [0, a] using erf functions
fprintf('Analytic method\n');
w = 2*a;
kappa = 0.01;
epsil = 0.0001;      % Small value for accuracy in polynomial degree calculation
k = sqrt(2)/kappa*sqrt(log(2/(pi*epsil^2)));
delta = (w+kappa)/2;
x = chebfun(@(x) x);
f_obj = chebfun({0, x/a, 0}, [-1 -a a 1]);  % Target function representation

% Degree of the polynomial for approximation
deg_list = [1001];  % Example degree for approximation

for j = 1:length(deg_list)
    deg = deg_list(j);
    fprintf('deg analytic = %d\n', deg);
    f_approx = x / (a) * (0.5 * (erf(k*(x+delta)) + erf(k*(-x+delta))));
    coef_full = chebcoeffs(fscale * f_approx, deg+1);
    poly_cheb = chebfun(coef_full, 'coeffs');
    max_cheb = max(abs(poly_cheb));
    if max_cheb > 1
        error('Inf norm of the function is > 1');
    end
end

% Plotting the analytic method results
figure(1)
clf
subplot(2,1,1);
hold on
plot(x, poly_cheb, 'ro-', 'LineWidth', 1.5);
plot(x, fscale*f_obj, 'b-', 'LineWidth', 2);
hold off
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 15);
legend({'Poly', 'Target'}, 'FontSize', 15);
title(sprintf('Analytic. Degree = %d', deg), 'FontSize', 15);
box on

subplot(2,1,2);
plot(x, abs(poly_cheb - fscale*f_obj), 'k-', 'LineWidth', 1.5);
xlim([0 a-0.01]);  % Limiting x-axis to [0, a]
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$|f_{\mathrm{poly}}(x) - f(x)|$', 'Interpreter', 'latex', 'FontSize', 15);
title('Error', 'FontSize', 15);

print(sprintf('uniform_singular_value_amplification_analytic_deg_%d.png', deg),'-dpng','-r500');


%% Convex Optimization Method
% Finds the optimal polynomial using convex optimization
fprintf('Convex optimization method\n');
deg_cvx_list = [51];  % Example degree list for convex optimization method
opts.npts = 1000;
opts.epsil = 1-fscale;
opts.fscale = fscale;
opts.intervals = [0, a];
opts.objnorm = Inf;
opts.isplot = false;

for j = 1:length(deg_cvx_list)
    deg = deg_cvx_list(j);
    fprintf('deg cvx = %d\n', deg);
    coef_full = cvx_poly_coef(targ, deg, opts);
    poly_cheb = chebfun(coef_full, 'coeffs');
    max_cheb = max(abs(poly_cheb));
    if max_cheb > 1
        error('Inf norm of the function is > 1');
    end
end

% Plotting the convex optimization method results
figure(2)
clf
subplot(2,1,1);
hold on
plot(x, poly_cheb, 'ro-', 'LineWidth', 1.5);
plot(x, fscale*f_obj, 'b-', 'LineWidth', 2);
hold off
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 15);
legend({'Poly', 'Target'}, 'FontSize', 15);
title(sprintf('Optimized. Degree = %d', deg), 'FontSize', 15);
box on

subplot(2,1,2);
plot(x, abs(poly_cheb - fscale*f_obj), 'k-', 'LineWidth', 1.5);
xlim([0 a-0.01]);  % Limiting x-axis to [0, a]
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$|f_{\mathrm{poly}}(x) - f(x)|$', 'Interpreter', 'latex', 'FontSize', 15);
title('Error', 'FontSize', 15);

print(sprintf('uniform_singular_value_amplification_convex_deg_%d.png', deg),'-dpng','-r500');

%% Find the QSP phase factors
% We use Newton's method for solving phase factors. The parameters of the 
% solver is initiated as follows.
opts.maxiter = 100;
opts.criteria = 1e-12;
opts.useReal = true;
opts.targetPre = true;
opts.method = 'Newton';
coef = coef_full(1+parity:2:end);
[phi_proc,out] = QSP_solver(coef,parity,opts);

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
print(sprintf('phase_uniform_singular_value_amplification_convex_deg_%d.png', deg),'-dpng','-r500');


%% References
% [1] Low, G. H., & Chuang, I. L. (2017). Hamiltonian Simulation by
% Uniform Spectral Amplification. http://arxiv.org/abs/1707.05391,
% [2] Gilyén, A., Su, Y., Low, G. H., & Wiebe, N. (2018). Quantum singular
% value transformation and beyond: exponential improvements for quantum
% matrix arithmetics. 1–67. http://arxiv.org/abs/1806.01838

