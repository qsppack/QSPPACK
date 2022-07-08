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
%                should be even/odd, only provide non-zero coefficients
%     parity --- Parity of polynomial P (0 -- even, 1 -- odd)
%       opts --- Options structure with fields
%                criteria: stop criteria 
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
% Author: Xiang Meng, Yulong Dong
% Version 1.0
% Last Update 06/2020
%
%--------------------------------------------------------------------------
% setup options for L-BFGS solver

if ~isfield(opts,'criteria');              opts.criteria = 1e-12; end
if ~isfield(opts,'useReal');               opts.useReal = true; end
if ~isfield(opts,'targetPre');             opts.targetPre = true;    end

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
[phi,obj_value,out] = QSP_LBFGS(obj,grad,delta,zeros(tot_len,1),opts);
% add pi/4 at the end
if (opts.targetPre == true)
    phi(end) = phi(end) + pi/4;
end
% construct full phase factors
if( parity == 0 )
  phi_proc = zeros(2*length(phi)-1,1);
  phi_proc(1:(length(phi)-1)) = phi(end:-1:2);
  phi_proc(length(phi):end)=phi;
else
  phi_proc = zeros(2*length(phi),1);
  phi_proc(1:length(phi)) = phi(end:-1:1);
  phi_proc(length(phi)+1:end) = phi;
end
time = toc;

%--------------------------------------------------------------------------
% output information

out.time = time;
out.value = obj_value;

end
