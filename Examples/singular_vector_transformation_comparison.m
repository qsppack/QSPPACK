%% Uniform Singular Value Amplification Comparison
% This script compares the performance of direct construction and convex optimization-based methods
% for uniform singular value amplification. It aims to approximate a linear function within a specified
% interval using these methods and evaluates their performance.

%% Initialization
fscale = 1/1.05;       % Scaling factor to ensure bounds
a = 0.2;            % Upper limit of the interval [0, a]
kappa = a * 2;
parity = 1;         % Ensures the function is odd
targ = @(x) sign(x);    % Target linear function scaled by 1/a

%% Analytic Method using Odd Polynomial Approximation
% Approximates the truncated linear function on [0, a] using erf functions
fprintf('Analytic method\n');
epsil = 0.0001;     % Small value for accuracy in polynomial degree calculation
k = sqrt(2)/kappa*sqrt(log(2/(pi*epsil^2)));
x = chebfun(@(x) x);
f_obj = sign(x);

% Degree of the polynomial for approximation
deg_list = [51];  % Example degree for approximation

for j = 1:length(deg_list)
    deg = deg_list(j);
    fprintf('deg analytic = %d\n', deg);
    coef_full = chebcoeffs(fscale * erf(k*x), deg+1);
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
xlim([a, 1]);  
ylim([0,2e-3])
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$|f_{\mathrm{poly}}(x) - f(x)|$', 'Interpreter', 'latex', 'FontSize', 15);
title('Error', 'FontSize', 15);

print(sprintf('singular_vector_amplification_analytic_deg_%d.png', deg),'-dpng','-r500');


%% Convex Optimization Method
% Finds the optimal polynomial using convex optimization
fprintf('Convex optimization method\n');
deg_cvx_list = [31];  % Example degree list for convex optimization method
opts.npts = 1000;
opts.epsil = 1-fscale;
opts.fscale = fscale;
opts.intervals = [a, 1];
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
xlim([a,1]); 
ylim([0,2e-3])
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$|f_{\mathrm{poly}}(x) - f(x)|$', 'Interpreter', 'latex', 'FontSize', 15);
title('Error', 'FontSize', 15);

print(sprintf('singular_vector_amplification_convex_deg_%d.png', deg),'-dpng','-r500');

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

xlist = linspace(a,1,500)';
func = @(x) ChebyCoef2Func(x, coef, parity, true);
targ_value = targ(xlist);
func_value = func(xlist);
QSP_value = QSPGetEntry(xlist, phi_proc, out);
err= norm(QSP_value-func_value,Inf);
disp('The residual error is');
disp(err);


%% References
% [1] Low, G. H., & Chuang, I. L. (2017). Hamiltonian Simulation by
% Uniform Spectral Amplification. http://arxiv.org/abs/1707.05391,
% [2] Gilyén, A., Su, Y., Low, G. H., & Wiebe, N. (2018). Quantum singular
% value transformation and beyond: exponential improvements for quantum
% matrix arithmetics. 1–67. http://arxiv.org/abs/1806.01838

