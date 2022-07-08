function [f,df] = F_Jacobian(phi, parity,opts)
%--------------------------------------------------------------------------
% Compute Pim and its Jacobian matrix values at sigle ponit cos(theta)
% P_im: the imagrinary part of the (1,1) element of the QSP unitary matrix.

% Input:
%       phi --- reduced phase factors 
%     theta --- vector
%    parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%         f --- F value at the point cos(theta)
%        df --- Jacobian matrix value at the point cos(theta)
%                  
%--------------------------------------------------------------------------
% setup options 
if ~isfield(opts,'useReal');              opts.useReal = true; end


%%--------------------------------------------------------------------------
% initial preparation
if (opts.useReal == true)
    f = @(x) QSPGetPimDeri_sym_real(phi, x, parity);
else
    f = @(x) QSPGetPimDeri_sym(phi, x, parity);
end
d = length(phi);
dd = 2 * d; 
theta = (0:d)'*pi/dd;
M = zeros(2*dd, d+1);

for n = 1:(d+1)
    M(n,:)=f(cos(theta(n)));
end

M(d+2:dd+1,:)=(-1)^parity * M(d:-1:1,:);
M(dd+2:end,:)= M(dd:-1:2,:);

M = fft(M); %FFT w.r.t. columns.
M = real(M(1:(dd+1),:));
M(2:end-1,:) = M(2:end-1,:)*2;
M = M/(2*dd);

f = M(parity+1:2:2*d,end);
df = M(parity+1:2:2*d,1:end-1);

end