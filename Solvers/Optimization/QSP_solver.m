function [phi_proc,out] = QSP_solver(coef,parity,opts)
%--------------------------------------------------------------------------
% Given coefficients of a polynomial P, yield corresponding phase factors
%
% The reference chose the first half of the phase factors as the 
% optimization variables, while in the code we used the second half of the 
% phase factors. These two formulations are equivalent. 
%
% To simplify the representation, a constant pi/4 is added to both sides of 
% the phase factors when evaluating the objective and the gradient. In the
% output, the FULL phase factors with pi/4 are given.
%
% Input:
%       coef --- Coefficients of polynomial P under Chevyshev basis, P
%                should be even/odd, only provide non-zero coefficients.
%                Coefficients should be ranked from low order term to high
%                order term.
%     parity --- Parity of polynomial P (0 -- even, 1 -- odd)
%       opts --- Options structure with fields
%                criteria: stop criteria  
%                useReal: use only real arithmetics if true
%                targetPre: want Pre to be target function if true
%                method: choose from 'LBFGS', 'FPI', or 'Newton'
%                typePhi: full or reduced phase factors
%
% Output:
%    phi_proc --- Solution of optimization problem, FULL phase factors
%         out --- Information of solving process
%
%--------------------------------------------------------------------------
%
% Reference: Yulong Dong, Xiang  Meng, K.Birgitta Whaley and Lin Lin
%            Efficient Phase Factor Evaluation in Quantum Signal Processing
%
% Version 2.0
% Last Update 08/2022
%
%--------------------------------------------------------------------------
% setup options for L-BFGS solver
if ~isfield(opts,'maxiter');               opts.maxiter = 5e4; end
if ~isfield(opts,'criteria');              opts.criteria = 1e-12; end
if ~isfield(opts,'useReal');               opts.useReal = true; end
if ~isfield(opts,'targetPre');             opts.targetPre = true;    end
if ~isfield(opts,'method');                opts.method = 'FPI'; end
if ~isfield(opts,'typePhi');               opts.typePhi = 'full'; end


if strcmp(opts.method,'LBFGS')
    
    %--------------------------------------------------------------------------
    % initial preparation
    
    tot_len = length(coef);
    delta = cos((1:2:(2*tot_len-1))*(pi/2/(2*tot_len)))';
    if (opts.targetPre == false)
        opts.target = @(x, opts) -ChebyCoef2Func(x, coef, parity, true);
    else
        opts.target = @(x, opts) ChebyCoef2Func(x, coef, parity, true);
    end
    opts.parity = parity;
    obj = @QSPObj_sym;
    if (opts.useReal == true)
        grad = @QSPGrad_sym_real;
    else
        grad = @QSPGrad_sym;
    end
    
    %--------------------------------------------------------------------------
    % solve by L-BFGS with selected initial point
    
    tic;
    [phi,err,iter] = QSP_LBFGS(obj,grad,delta,zeros(tot_len,1),opts);
    % convert phi to reduced phase factors
    if( parity == 0 )
        phi(1) = phi(1)/2;
    end
    runtime = toc;

elseif strcmp(opts.method,'FPI')
    [phi, err, iter, runtime] = QSP_CM(coef, parity, opts);
    
elseif strcmp(opts.method, 'Newton')
    [phi, err, iter, runtime] = QSP_Newton(coef, parity, opts);
else
    fprintf("Assigned method doesn't exist. Please choose method from 'LBFGS', 'FPI' or 'Newton'.\n");
end

%--------------------------------------------------------------------------
% output information

out.iter = iter;
out.time = runtime;
out.value = err;
out.parity = parity;
out.targetPre = opts.targetPre;

if strcmp(opts.typePhi,'full')
    phi_proc = rdc_phase_factor_to_full(phi, parity, opts.targetPre);
    out.typePhi = 'full';
else
    phi_proc = phi;
    out.typePhi = 'reduced';
end

end
