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

d = length(phi);
expphi = exp(1j*phi);
Wx = [x 1j*sqrt(1-x^2); 1j*sqrt(1-x^2) x];
right = zeros(2,d);
left = zeros(d,2);
y = zeros(1,d+1);
right(:,end) = [expphi(end);0];

for k = d-1:-1:1
    right(:,k) = [expphi(k); conj(expphi(k))] .* Wx * right(:,k+1);
end
left(1,:) = transpose(right(:,1));

if (parity==1)
    left(1,:) = left(1,:) * Wx;
end
for k = 1:d-1
    left(k+1,:) =  left(k,:) .* [expphi(k) conj(expphi(k))] * Wx;
end

for k = 1:d
    y(k) = 2 * left(k,:) .* [1j -1j] * right(:, k);
end

y(d+1) = left(d,:) * right(:,d);
y = imag(y);
end