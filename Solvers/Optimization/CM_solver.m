function [phi_full, step_err, step_count, runtime] = CM_solver(coef, parity, opts)
%--------------------------------------------------------------------------
% Contraction mapping solver for finding phase factors such that the real 
% part of the (1,1) element of the QSP unitary matrix gives desire Chebyshev
% expansion.

% Input:
%           coef --- Chebyshev coefficients
%         parity --- Parity of phi (0 -- even, 1 -- odd)
%           opts --- Options structure with fields
%                   maxiter: maximal iteration number
%                   criteria: stop criteria 
%                   targetPre: Pre to be target function 
%                   useReal: use real matrix mulplitcation to get QSP entry
%
% Output:
%       phi_full --- full set of phase factors
%       step_err --- error (1 norm)
%     step_count --- iteration number
%        runtime --- time used 
%                  
%--------------------------------------------------------------------------
% setup options for CM solver
if ~isfield(opts,'maxiter');              opts.maxiter = 1e5; end
if ~isfield(opts,'criteria');             opts.criteria = 1e-12; end
if ~isfield(opts,'targetPre');            opts.targetPre = true;    end
if ~isfield(opts,'useReal');              opts.useReal = true; end


tic

%%--------------------------------------------------------------------------
% initial preparation
if (opts.targetPre == true)   
    coef = - coef; % inverse is necessary
end
phi = zeros(size(coef));
Fval = F(phi,parity,opts);
step_err = norm(coef - Fval, 1);
step_count = 0;


%--------------------------------------------------------------------------
% solve by contraction mapping algorithm

while step_err > opts.criteria && step_count < opts.maxiter
    phi = phi + coef/2 - Fval/2;
    Fval = F(phi,parity,opts);
    step_err = norm(coef - Fval, 1);
    step_count = step_count + 1;
end

phi_full = rdc_phase_factor_to_full(phi, parity,opts.targetPre);
runtime = toc;
end