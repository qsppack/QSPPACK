function [phi, err, runtime] = QSP_Weiss_INFFT(coef, parity, opts)
%--------------------------------------------------------------------------
% A fast solver for finding phase factors such that the real part of the 
% (1,1) element of the QSP unitary matrix gives desire Chebyshev expansion.
% The solver first uses the Weiss algorithm to compute the complimentary polynomial,
% and then applies the inverse nonlinear fast Fourier transform (INFFT) to calculate
% the phase factors.

% Reference:
% Ni, H., Sarkar, R., Ying, L., & Lin, L. (2025). Inverse nonlinear fast Fourier
% transform on SU (2) with applications to quantum signal processing. arXiv preprint
% arXiv:2505.12615.

% Input:
%         coef --- Chebyshev coefficients
%       parity --- Parity of phi (0 -- even, 1 -- odd)
%          opts --- Options structure with fields
%                   Weiss_thd: the target accuracy for Weiss algorithm 
%                   eta: An estimation of the lower bound of 1 - max_{|x|<=1}(|f(x)|). A wrong value (e.g. the default value 0.5) won't affect the result, but may slightly affect the speed.
%                   targetPre: Pre to be target function 
%
% Output:
%       phi --- reduced phase factors
%       err --- error (2 norm)
%   runtime --- time used 
%  
% ----------------------------------------------------------------------
% Author:    Hongkang Ni   update 06/2025
%
%--------------------------------------------------------------------------

if ~isfield(opts,'eta');                opts.eta = 0.5; end
if ~isfield(opts,'Weiss_thd');          opts.Weiss_thd = 1e-10; end
if ~isfield(opts,'targetPre');            opts.targetPre = true;    end

tic
% Transpose if coef is a row vector
if ~isvector(coef)
    error('invalid input')
end
flag = 0;
if isrow(coef)
    coef = coef.';
    flag = -1;
end

if (opts.targetPre == true)   
    coef = - coef; % inverse is necessary
end


%% Weiss algorithm
if parity == 0     
    coef(1) = coef(1)*2;   
end
bc = coef/2;
d = length(bc) - 1;
N = ceil(d/opts.eta);
N = max(N,2*d+2);
thd = 1;

while thd > opts.Weiss_thd
    ext_bc = [bc;zeros(N-2*d-1-parity,1);bc(end:-1:2-parity)];
    bz = ifft(ext_bc)*N;
    if parity == 1
        bz = bz.*exp(1i*pi/N*(0:N-1)');  
    end
    logsqrt_b = log(sqrt(1-abs(bz).^2));
    clear ext_bc bz;
    r = fft(logsqrt_b)/N;
    r(2:N/2)=0;
    r(N/2+1:N) = 2*r(N/2+1:N);
    G = ifft(r)*N;
    clear logsqrt_b r;
    AN = fft(conj(exp(G)))/N;
    thd = norm(AN(floor(N/4):ceil(N*3/4)))^2/norm(AN)^2;
    N = N * 3;
end


%% inverse nonlinear fast Fourier transform (INFFT)
ac = real(AN(1:d+1)); % coefficient of a^*(z)
[F2, xi1, eta1] = INFFT(ac,flip(bc));
phi = atan(F2(end:-1:1));
runtime = toc;

%% error calculation
xi2 = circshift(flip(xi1),-1);
eta2 = eta1;
if parity == 0
    gamma0 = F2(end);
    xi_temp = (xi2 - gamma0*eta2)/sqrt(1+abs(gamma0)^2);
    eta2 = (gamma0*xi2 + eta2)/sqrt(1+abs(gamma0)^2);
    xi2 = xi_temp(2:end);
    eta2 = eta2(1:end-1);
end
xin = convfft(xi1, eta2) + convfft(flip(eta1), xi2);
xin = xin(1:d+1);
err = norm(xin*2-flip(coef));


if flag == -1
    phi = phi.';
end

