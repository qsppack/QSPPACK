function Pim = QSPGetPim_sym(phi, y, parity)
%--------------------------------------------------------------------------
% Get the imaginary value of polynomial P in the QSP unitary matrix based 
% on given reduced phase vector and point y \in [-1, 1] 
% y can be a list of point, and in this case P will also be a list.

% Input:
%       phi --- reduced phase factors   
%    parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%       Pim --- The imaginary part of the (1,1) element of the QSP unitary 
%               matrix at point y.
%--------------------------------------------------------------------------
% Author:     Jiasu Wang   update 07/2022
% 
%--------------------------------------------------------------------------

Pim = y;
phi = phi(end:-1:1);
expphi = exp(1j*phi);
for n = 1:length(y)
    x = y(n);
    Wx = [x 1j*sqrt(1-x^2); 1j*sqrt(1-x^2) x];
    
    ret = [expphi(1) 0];
    for k = 2:length(expphi)
        ret = ret*Wx*[expphi(k) 0;0 conj(expphi(k))];
    end
    if(parity==1)
        P = ret * Wx * transpose(ret); 
    else
        P = ret * transpose(ret); 
    end
    Pim(n) = imag(P);
end
end
