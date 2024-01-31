%% Generate the phase factors for a random Chebyshev polynomial. 
%  
%  The purpose of this example is to verify that the modified phase factor
%  generated using the subroutine modify_phase_factor_circuit correctly
%  encodes the polynomial when evaluated using
%  QSPGetUnitary_modified_phase.
%
%  This needs to test for all cases d=4k+1,4k+2,4k+3,4k+4, as well as the
%  negation of the phase factors to verify the conjugation relation.

%%
% (example/chebyshev_random_modified_phase.m)

deg_list = [9,10,11,12];

for deg = deg_list
    disp('deg = ')
    disp(deg);
    parity = mod(deg,2);
    coef_targ = randn(deg+1,1);
    if parity == 0
      % Eliminate odd terms
      coef_targ(2:2:end) = 0;
    else
      % Eliminate even terms
      coef_targ(1:2:end) = 0;
    end
    targ = chebfun(coef_targ, 'coeffs');
    max_targ = max(abs(targ));
    % normalize and regenerate the Chebyshev polynomial
    coef_targ = coef_targ / max_targ * 0.9;
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
    % Convert to the phase factor used in the circuit
    phi_mod=modify_phase_factor_circuit(phi_proc);
    
    % We do the following test to demonstrate that the obtained phase factors 
    % satisfy expectation.
    xlist = linspace(-1, 1, 100)';
    targ_value = targ(xlist);
    QSP_value_pos = zeros(size(xlist));
    QSP_value_neg = zeros(size(xlist));    
    for j = 1 : length(xlist)
        x = xlist(j);
        QSP_value_pos(j) = QSPGetUnitary_modified_phase(phi_mod, x);
        QSP_value_neg(j) = QSPGetUnitary_modified_phase(-phi_mod, x);
    end
    err1= norm(real(QSP_value_pos)-targ_value,1)/length(xlist);
    disp('The error of real(QSP_pos) is');
    disp(err1);
    err2= norm(0.5*(QSP_value_pos+QSP_value_neg)-targ_value,1)/length(xlist);
    disp('The error of (QSP_pos+QSP_neg)/2 is');
    disp(err2);
    
    pause
end
