%% Generate the phase factors for a Chebyshev polynomial
% Phase factor for a simple 3rd order Chebyshev polynomial
%
% $f(x) = 0.2 T_1(x) + 0.4 T_3(x)$

%%
% (example/chebyshev_ex3.m)

deg = 3;
parity = mod(deg,2);
coef_targ = zeros(deg+1,1);
coef_targ(2)=0.2;
coef_targ(4)=0.4;
targ = chebfun(coef_targ, 'coeffs');

%%
% We only need its Chebyshev coefficients with respect to its parity
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

%% 
% Let us check it explicitly (i.e., we look into the implementation of
% QSPGetUnitary)

x = rand(1);
fprintf('at a random x = %g\n', x);
fprintf('QSP value = %g\n', QSPGetUnitary(phi_proc,x));
fprintf('Target    = %g\n', targ(x));

function [targ] = QSPGetUnitary(phase, x)

    Wx = [x, 1j*sqrt(1-x^2); 1j*sqrt(1-x^2), x];
    expphi = exp(1j*phase);

    ret = [expphi(1), 0; 0, conj(expphi(1))];

    for k = 2:numel(expphi)
        temp = [expphi(k), 0; 0, conj(expphi(k))];
        ret = ret * Wx * temp;
    end

    targ = real(ret(1,1));

end




