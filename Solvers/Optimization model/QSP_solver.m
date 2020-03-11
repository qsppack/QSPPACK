function [phi,out] = QSP_solver(coef,parity,opts)
%--------------------------------------------------------------------------
% Given coefficients of a polynomial P, yield corresponding phase factors
%
% The reference chose the first half of the phase factors as the 
% optimization variables, while in the code we used the second half of the 
% phase factors. These two formulations are equivalent. 
%
% In addition, a constant pi/4 is added to both sides of the phase factors 
% to simplify the representation.
%
% Input:
%       coef --- Coefficients of polynomial P under Chevyshev basis, P
%                should be even/odd, only provide non-zero coefficients
%     parity --- Parity of polynomial P (0 -- even, 1 -- odd)
%       opts --- Options structure with fields
%                criteria: stop criteria 
%
% Output:
%         phi --- Solution of robust optimization problem
%         out --- Information of solving process
%
%--------------------------------------------------------------------------
%
% Reference: Yulong Dong, Xiang  Meng, K.Birgitta Whaley and Lin Lin
%            Efficient Phase Factor Evaluation in Quantum Signal Processing
%
% Author: Xiang Meng, Yulong Dong
% Version 1.0
% Last Update 2020.2
%
%--------------------------------------------------------------------------
% setup options for L-BFGS solver

if ~isfield(opts,'criteria');              opts.criteria = 1e-12; end

%--------------------------------------------------------------------------
% initial preparation

tot_len = length(coef);
delta = cos((1:2:(2*tot_len-1))*(pi/2/(2*tot_len)))';
opts.target = @(x, opts) ChebyCoef2Func(x, coef, 1-parity, true);
opts.parity = parity;
obj = @QSPObj_sym;
grad = @QSPGrad_sym;

%--------------------------------------------------------------------------
% solve by L-BFGS with selected initial point
% we've addeded \pi/4 on both sides in function obj and grad 

tic;
[phi,obj_value,out] = QSP_LBFGS(obj,grad,delta,zeros(tot_len,1),opts);
time = toc;

%--------------------------------------------------------------------------
% output information

out.time = time;
out.value = obj_value;

end
