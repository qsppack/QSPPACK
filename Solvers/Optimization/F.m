function coe = F(phi, parity, opts)
%--------------------------------------------------------------------------
% Compute the Chebyshev coefficients of P_im 
% P_im: the imagrinary part of the (1,1) element of the QSP unitary matrix.

% Input:
%          phi --- reduced phase factors 
%       parity --- Parity of phi (0 -- even, 1 -- odd)
%
% Output:
%          coe --- Chebyshev coefficients of P_im w.r.t. 
%                  T_(2k) -- even, or T_(2k-1) -- odd
% ----------------------------------------------------------------------
% Author:     Hongkang Ni   update 04/2022
%             Jiasu Wang    update 07/2022
%                  
%--------------------------------------------------------------------------
% setup options for CM solver
if ~isfield(opts,'useReal');              opts.useReal = true; end


%%--------------------------------------------------------------------------
% initial preparation
d = length(phi);
dd = 2*d;
theta = (0:d)'*pi/dd;
M = zeros(2*dd, 1);

if (opts.useReal == true)
    f = @(x) QSPGetPim_sym_real(phi, x, parity);
else
    f = @(x) QSPGetPim_sym(phi, x, parity);
end

%-------------------------------------------------------------------------
% start chebyshev coefficients evaluation

M(1:d+1) = f(cos(theta));
M(d+2:dd+1)=(-1)^parity * M(d:-1:1);
M(dd+2:end)= M(dd:-1:2);
M = fft(M); %FFT w.r.t. columns.
M = real(M);
M = M/(2*dd);
M(2:end-1) = M(2:end-1)*2;
coe = M(parity+1:2:2*d);

%--------------------------------------------------------------------------
end



