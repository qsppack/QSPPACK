function [gamma_list, xin, etan] = INFFT(ac, bc)
%--------------------------------------------------------------------------
% Compute the inverse NLFT of $(a(z),b(z))$ given the coefficients of
% $a^*(z)$ and $b(z)$ in $O(n log^2(n))$ time. See details in the paper
% Ni, H., Sarkar, R., Ying, L., & Lin, L. (2025). Inverse nonlinear fast Fourier
% transform on SU (2) with applications to quantum signal processing. arXiv preprint
% arXiv:2505.12615.
%
% Input:
%       ac --- (column vector) coefficient list of $a^*(z)$
%       bc --- (column vector) coefficient list of $b(z)$
%
% Output:
%         gamma_list --- the NLFT coefficients 
%         xin, etan --- helper variables for recursion, not used in the final output
%  
% ----------------------------------------------------------------------
% Author:    Hongkang Ni   update 06/2025
%
%--------------------------------------------------------------------------

n = length(ac);
gamma_list = zeros(n,1);
if n == 1
    Gi = [ac,bc];
    gamma = Gi(1,2)/Gi(1,1);
    gamma_list(1) = gamma;
    xin = [gamma_list(1);0]/sqrt(1+abs(gamma_list(1))^2);
    etan = [1;0]/sqrt(1+abs(gamma_list(1))^2);
else
    m = ceil(n/2);
    [gamma_first_half, xi1, eta1] = INFFT(ac(1:m), bc(1:m));
    gamma_list(1:m) = gamma_first_half;
    
    am = convfft(eta1(end:-1:1),ac,true) + convfft(xi1(end:-1:1),bc,true);
    bm = -convfft(xi1,ac,true) + convfft(eta1,bc,true);
    
    [gamma_last_half, xi2, eta2] = INFFT(am, bm);
    gamma_list(m+1:n) = gamma_last_half;
    
    xin = convfft(xi1, eta2) + convfft(eta1(end:-1:1), xi2);
    etan = convfft(eta1, eta2) - convfft(xi1(end:-1:1), xi2);
end
    
    
    