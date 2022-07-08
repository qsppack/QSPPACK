function Pim = QSPGetPim_sym_real(phi, y, parity)
%--------------------------------------------------------------------------
% Get the imaginary value of polynomial P in the QSP unitary matrix based 
% on given reduced phase vector and point y \in [-1, 1] 
% using the real matrix representation of P_im
% y can be a list of point, and in this case P will also be a list.

% Input:
%       phi --- reduced phase factors   
%    parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%       Pim --- The imaginary part of the (1,1) element of the QSP unitary 
%               matrix at point y.
%--------------------------------------------------------------------------

Pim = y;
n = length(phi);
for m = 1:length(y)
    theta = acos(y(m));
    B = [cos(2*theta), 0, -sin(2*theta);
        0, 1, 0;
        sin(2*theta), 0, cos(2*theta)];
    
    if parity == 0
        R = [1;0;0];
    else
        R = [cos(theta); 0; sin(theta)];
    end
    for k = 2:n
        R = B*[cos(2*phi(k-1)), -sin(2*phi(k-1)), 0; sin(2*phi(k-1)), cos(2*phi(k-1)), 0; 0, 0, 1]*R;
    end
    
    Pim(m) =[sin(2*phi(n)), cos(2*phi(n)), 0]*R;
end
end