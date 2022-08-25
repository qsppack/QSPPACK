function qspmat = QSPGetUnit_sym(phi, x, parity)
%--------------------------------------------------------------------------
% Get the QSP unitary matrix based on given phase vector and point 
% x \in [-1, 1]
%
% Input:
%       phi --- The last half of phase factors (If parity == 1, then phi is
%               the reduced phase factor. If parity == 0, then phi(1) is
%               different from the reduced phase factors by a factor of 2.)
%         x --- Point to be evaluated
%    parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%     qspmat--- The QSP unitary matrix
%
%--------------------------------------------------------------------------
% Author:    Xiang Meng, Yulong Dong  update 02/2020
%
%--------------------------------------------------------------------------
% compute the QSP matrix

Wx = [x 1j*sqrt(1-x^2); 1j*sqrt(1-x^2) x];
gate = [exp(1j*pi/4) 0;0 conj(exp(1j*pi/4))];
expphi = exp(1j*phi);

if(parity==1)
    ret = [expphi(1) 0;0 conj(expphi(1))];
    for k = 2:length(expphi)
        ret = ret*Wx*[expphi(k) 0;0 conj(expphi(k))];
    end
    ret = ret*gate;
    qspmat = transpose(ret)*Wx*ret;
else
    ret = eye(2);
    for k = 2:length(expphi)
        ret = ret*Wx.*[expphi(k) conj(expphi(k))];
    end
    ret = ret*gate;
    qspmat = transpose(ret)*[expphi(1) 0;0 conj(expphi(1))]*ret; 
end

%--------------------------------------------------------------------------

end
