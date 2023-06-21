function [grad, obj] = QSPGrad_sym_real(phi,delta,opts)
%--------------------------------------------------------------------------
% Evalute the gradient and objective of QSP function, provided that 
% phi is symmetric
%
% Input:
%        phi --- Variables
%      delta --- Samples
%       opts --- Options structure with fields
%         target: target function
%         parity: parity of phi (0 -- even, 1 -- odd)
%
% Output:
%        grad --- Gradient of obj function
%         obj --- Objective function value
%
% ----------------------------------------------------------------------
% Author:    Jiasu Wang  update 07/2022
%
%--------------------------------------------------------------------------
% initial computation

m = length(delta);
d = length(phi);
obj = zeros(m,1);
grad = zeros(m,d);
targetx = opts.target;
parity = opts.parity;

% convert the phase factor used in LBFGS solver to reduced phase factors
if (parity==0)
    phi(1)=phi(1)/2; 
end

%-------------------------------------------------------------------------
% start gradient evaluation

for i=1:m
    x = delta(i); 
    y = QSPGetPimDeri_sym_real(phi, x, parity);
    if (parity==0)
        y(1)=y(1)/2;
    end
    y = -y;% flip the sign
    gap = y(end)-targetx(x);
    obj(i) = 0.5 * gap^2;
    grad(i,:)=y(1:end-1)*gap;
end

%--------------------------------------------------------------------------

end   