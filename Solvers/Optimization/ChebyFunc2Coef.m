function c = ChebyFunc2Coef(func, maxorder)
%--------------------------------------------------------------------------
% Evaluate Chebyshev coefficients of a real polynomial of degree at most maxorder            
%--------------------------------------------------------------------------
%
% Input:
%       func, maxorder
% Output:
%       c --- Chebyshev coefficients up to maxorder
%
% ----------------------------------------------------------------------

    M = maxorder;       
    theta = zeros(2*M,1);
    for i=1:2*M
        theta(i) = (i-1)*pi/M;
    end
    f = func(cos(theta));
    c = fft(f);
    c = real(c);
    c = c(1:M+1);
    c(2:end-1) = c(2:end-1)*2;
    c = c/(2*M);
end