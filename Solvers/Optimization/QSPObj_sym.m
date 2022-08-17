function [obj] = QSPObj_sym(phi, delta, opts)
%--------------------------------------------------------------------------
% Evalute the objective of QSP function, provided that phi is symmetric.
%
% Input:
%        phi --- Variables (If parity == 1, then phi is the reduced phase
%                factor. If parity == 0, then phi(1) is different from the
%                reduced phase factors by a factor of 2.)
%      delta --- Samples
%       opts --- Options structure with fields
%         target: target function
%         parity: parity of phi (0 -- even, 1 -- odd)
%
% Output:
%         obj --- Objective function value
%
%--------------------------------------------------------------------------
%
% Author: Xiang Meng, Yulong Dong
% Version 1.0 
% Last revision 02/2020
%
%--------------------------------------------------------------------------
% compute the objective

m = length(delta);
obj = zeros(m,1);
for i = 1:m
    qspmat = QSPGetUnit_sym(phi, delta(i), opts.parity);
    obj(i) = 0.5*(real(qspmat(1,1))-opts.target(delta(i)))^2;
end

%--------------------------------------------------------------------------

end
