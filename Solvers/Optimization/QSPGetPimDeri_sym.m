function y = QSPGetPimDeri_sym(phi, x, parity)
%--------------------------------------------------------------------------
% Compute Pim and its Jacobian matrix values at sigle ponit x
% P_im: the imagrinary part of the (1,1) element of the QSP unitary
% matrix.
% Note: theta MUST be a number

% Input:
%          phi --- reduced phase factors 
%        theta --- vector
%       parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%            y --- Pim and its Jacobian matrix value at the point x
%                  
%--------------------------------------------------------------------------
% Author:     Jiasu Wang   update 07/2022
% 
%--------------------------------------------------------------------------

d = length(phi);
expphi = exp(1j*phi);
Wx = [x 1j*sqrt(1-x^2); 1j*sqrt(1-x^2) x];
right = zeros(2,d);
y = zeros(1,d+1);
right(:,end) = [expphi(end);0];

for k = d-1:-1:1
    right(:,k) = [expphi(k); conj(expphi(k))] .* Wx * right(:,k+1);
end
left = transpose(right(:,1));
right = [1j; -1j] .* right;

if (parity==1)
    left = left * Wx;
end

for k = 1:d-1
    y(k) = 2 * left * right(:, k);
    left =  left .* [expphi(k) conj(expphi(k))] * Wx;
end
y(d) = 2 * left * right(:, d);
y(d+1) = left * [expphi(end); 0];
y = imag(y);
end