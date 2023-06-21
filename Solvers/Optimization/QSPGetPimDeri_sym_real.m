function y = QSPGetPimDeri_sym_real(phi, x, parity)
%--------------------------------------------------------------------------
% Compute Pim and its Jacobian matrix values at single ponit x
% P_im: the imagrinary part of the (1,1) element of the QSP unitary matrix.
% Using the real matrix representation of Pim
% Note: theta MUST be a number
%
% Input:
%          phi --- reduced phase factors 
%        theta --- vector
%       parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%            y --- Pim and its Jacobian matrix value at the point x
%                  
%--------------------------------------------------------------------------
% Author:      Hongkang Ni  update 04/2022
%              Jiasu Wang   update 07/2022
% 
%--------------------------------------------------------------------------

n = length(phi);
theta = acos(x);
B = [cos(2*theta), 0, -sin(2*theta);
    0, 1, 0;
    sin(2*theta), 0, cos(2*theta)];
L = zeros(n,3);     
L(n,:) = [0,1,0];
for k = n-1:-1:1
    L(k,:) = L(k+1,:)*[cos(2*phi(k+1)), -sin(2*phi(k+1)), 0; sin(2*phi(k+1)), cos(2*phi(k+1)), 0; 0, 0, 1]*B;
end        
R = zeros(3,n);    
if parity == 0
    R(:,1) = [1;0;0];
else
    R(:,1) = [cos(theta); 0; sin(theta)];
end
for k = 2:n
    R(:,k) = B*([cos(2*phi(k-1)), -sin(2*phi(k-1)), 0; sin(2*phi(k-1)), cos(2*phi(k-1)), 0; 0, 0, 1]*R(:,k-1));
end

y = zeros(1,n);
for k = 1:n
    y(k) = 2*L(k,:)*[-sin(2*phi(k)), -cos(2*phi(k)), 0; cos(2*phi(k)), -sin(2*phi(k)), 0; 0, 0, 0]*R(:,k);
end
y(n+1) = L(n,:)*[cos(2*phi(n)), -sin(2*phi(n)), 0; sin(2*phi(n)), cos(2*phi(n)), 0; 0, 0, 1]*R(:,n);

end