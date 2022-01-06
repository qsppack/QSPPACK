function [hess, grad, obj] = QSPHess_sym(phi,delta,opts)
%--------------------------------------------------------------------------
% Evalute the hessian, gradient and objective of QSP function, provided
% that phi is symmetric
%
% Input:
%        phi --- Variables
%      delta --- Samples
%       opts --- Options structure with fields
%         target: target function
%         parity: parity of phi (0 -- even, 1 -- odd)
%
% Output:
%     hessian --- Hessian of obj function
%        grad --- Gradient of obj function
%         obj --- Objective function value
%
%--------------------------------------------------------------------------
%
% Author: Xiang Meng, Yulong Dong
% Version 1.0
% Last revision 02/2020
%
%--------------------------------------------------------------------------
% initial computation

m = length(delta);
d = length(phi);
obj = zeros(m,1);
grad = zeros(m,d);
hess = zeros(d,d,m);
sigmaz = [1 0;0 -1];
tmp_sig = 1i*sigmaz;
gate = [exp(1j*pi/4) 0;0 conj(exp(1j*pi/4))];
exptheta = exp(1j*phi);
targetx = opts.target;
parity = opts.parity;

%--------------------------------------------------------------------------
% start hessian evaluation

for i=1:m
    x = delta(i);
    Wx = [x, 1j*sqrt(1-x^2); 1j*sqrt(1-x^2), x];
    tmp_save = zeros(2,2,d,d+1);
    tmp_lefsave = zeros(2,2,d+1);
    for j=1:d
        tmp_save(:,:,j,j) = eye(2);
    end
    thetawx = zeros(2,2,d);
    for j=1:d
        thetawx(:,:,j) = [exptheta(j) 0;0 conj(exptheta(j))]*Wx;
    end
    for j=1:d
        for k=j+1:d+1
            if(k<d+1)
                tmp_save(:,:,j,k) = tmp_save(:,:,j,k-1)*thetawx(:,:,k-1);
            else
                tmp_save(:,:,j,k) = tmp_save(:,:,j,k-1)*[exptheta(d) 0;0 conj(exptheta(d))]*gate;
            end
        end
    end
    if(parity==1)
        qspmat = transpose(tmp_save(:,:,1,d+1))*Wx*tmp_save(:,:,1,d+1);
        gap = real(qspmat(1,1))-targetx(x);
        leftmat = transpose(tmp_save(:,:,1,d+1))*Wx;
        for j=1:d+1
            tmp_lefsave(:,:,j) = leftmat*tmp_save(:,:,1,j);
        end
        grad_tmp = zeros(2,2,d);
        grad_leftmp = zeros(2,2,d);
        grad_Wtmp = zeros(2,2,d);
        for j=1:d
            grad_tmp(:,:,j) = tmp_save(:,:,1,j)*tmp_sig*tmp_save(:,:,j,d+1);
            grad_leftmp(:,:,j) = 2*leftmat*grad_tmp(:,:,j);
            grad(i,j) = real(grad_leftmp(1,1,j))*gap;
        end
        for j=1:d
            grad_Wtmp(:,:,j) = Wx*grad_tmp(:,:,j);
        end
        for j=1:d
            for k=j:d 
                hesst = tmp_lefsave(:,:,j).*[1 -1]*tmp_save(:,:,j,k).*[-1 1]*tmp_save(:,:,k,d+1);
                hesstmp = 2*real(hesst(1,1))*gap;
                hesst2 = transpose(grad_tmp(:,:,j))*grad_Wtmp(:,:,k);
                hesstmp2 = 2*real(hesst2(1,1))*gap;
                hess(j,k,i) = real(grad_leftmp(1,1,j))*real(grad_leftmp(1,1,k))+hesstmp+hesstmp2;
                hess(k,j,i) = hess(j,k,i);
            end
        end
        obj(i) = 0.5*(real(qspmat(1,1))-targetx(x))^2;
    else
        qspmat = transpose(tmp_save(:,:,2,d+1))*Wx*tmp_save(:,:,1,d+1);
        gap = real(qspmat(1,1))-targetx(x);
        leftmat = transpose(tmp_save(:,:,2,d+1))*Wx;
        for j=1:d+1
            tmp_lefsave(:,:,j) = leftmat*tmp_save(:,:,1,j);
        end
        grad_tmp = zeros(2,2,d);
        grad_tmp2 = zeros(2,2,d);
        grad_leftmp = zeros(2,2,d);
        grad_Wtmp = zeros(2,2,d);
        for j=1:d
            if(j~=1)
                grad_tmp2(:,:,j) = tmp_save(:,:,2,j)*tmp_sig*tmp_save(:,:,j,d+1);
            end
            grad_tmp(:,:,j) = tmp_save(:,:,1,j)*tmp_sig*tmp_save(:,:,j,d+1);
            grad_leftmp(:,:,j) = 2*leftmat*grad_tmp(:,:,j);
            if j==1; grad_leftmp(:,:,j) = grad_leftmp(:,:,j)/2; end
            grad(i,j) = real(grad_leftmp(1,1,j))*gap;
        end
        for j=1:d
            grad_Wtmp(:,:,j) = Wx*grad_tmp(:,:,j);
        end
        for k=1:d
            hesst = leftmat.*[1 -1]*tmp_save(:,:,1,k).*[-1 1]*tmp_save(:,:,k,d+1);
            hesstmp = real(hesst(1,1))*gap;
            if(k==1)
                hesst2 = zeros(2,2);
            else
                hesst2 = transpose(tmp_save(:,:,1,d+1))*tmp_sig*Wx*grad_tmp2(:,:,k);
            end
            hesstmp2 = real(hesst2(1,1))*gap;
            hess(1,k,i) = real(grad_leftmp(1,1,1))*real(grad_leftmp(1,1,k))+hesstmp+hesstmp2;
            hess(k,1,i) = hess(1,k,i);
        end
        for j=2:d
            for k=j:d 
                hesst = tmp_lefsave(:,:,j).*[1 -1]*tmp_save(:,:,j,k).*[-1 1]*tmp_save(:,:,k,d+1);
                hesstmp = 2*real(hesst(1,1))*gap;
                hesst2 = transpose(grad_tmp2(:,:,j))*grad_Wtmp(:,:,k);
                hesstmp2 = 2*real(hesst2(1,1))*gap;
                hess(j,k,i) = real(grad_leftmp(1,1,j))*real(grad_leftmp(1,1,k))+hesstmp+hesstmp2;
                hess(k,j,i) = hess(j,k,i);
            end
        end
        obj(i) = 0.5*(real(qspmat(1,1))-targetx(x))^2;
    end
end

%--------------------------------------------------------------------------

end
