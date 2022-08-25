function [grad, obj] = QSPGrad_sym(phi,delta,opts)
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
%--------------------------------------------------------------------------
% Author:    Xiang Meng, Yulong Dong  update 02/2020
%
%--------------------------------------------------------------------------
% initial computation

m = length(delta);
d = length(phi);
obj = zeros(m,1);
grad = zeros(m,d);
gate = [exp(1j*pi/4) 0;0 conj(exp(1j*pi/4))];
exptheta = exp(1j*phi);
targetx = opts.target;
parity = opts.parity;

%-------------------------------------------------------------------------
% start gradient evaluation

for i=1:m
    x = delta(i);
    Wx = [x 1j*sqrt(1-x^2);1j*sqrt(1-x^2) x];
    tmp_save1 = zeros(2,2,d);
    tmp_save2 = zeros(2,2,d);
    tmp_save1(:,:,1) = eye(2);
    tmp_save2(:,:,1) = [exptheta(d) 0;0 conj(exptheta(d))]*gate;
    for j=2:d
        tmp_save1(:,:,j) = tmp_save1(:,:,j-1).*[exptheta(j-1) conj(exptheta(j-1))]*Wx;
        tmp_save2(:,:,j) = [exptheta(d-j+1);conj(exptheta(d-j+1))].*Wx*tmp_save2(:,:,j-1);
    end
    if(parity==1)
        qspmat = transpose(tmp_save2(:,:,d))*Wx*tmp_save2(:,:,d);
        gap = real(qspmat(1,1))-targetx(x);
        leftmat = transpose(tmp_save2(:,:,d))*Wx;
        for j=1:d
            grad_tmp = leftmat*tmp_save1(:,:,j).*[1j -1j]*tmp_save2(:,:,d+1-j);
            grad(i,j) = 2*real(grad_tmp(1,1))*gap;
        end
        obj(i) = 0.5*(real(qspmat(1,1))-targetx(x))^2;
    else
        qspmat = transpose(tmp_save2(:,:,d-1))*Wx*tmp_save2(:,:,d);
        gap = real(qspmat(1,1))-targetx(x);
        leftmat = transpose(tmp_save2(:,:,d-1))*Wx;
        for j=1:d
            grad_tmp = leftmat*tmp_save1(:,:,j).*[1j -1j]*tmp_save2(:,:,d+1-j);
            grad(i,j) = 2*real(grad_tmp(1,1))*gap;
        end
        grad(i,1) = grad(i,1)/2;
        obj(i) = 0.5*(real(qspmat(1,1))-targetx(x))^2;
    end
end

%--------------------------------------------------------------------------

end   
