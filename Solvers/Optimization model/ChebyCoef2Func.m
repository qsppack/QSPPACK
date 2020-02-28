% Evaluate function based on Chebyshev expansion coefficients
%
% ----------------------------------------------------------------------
%
% Input:
%       x, coef
%       parity  -- true for even
%  partialcoef  -- true: only include even/odd coefficiennts
% Output:
%       ret     -- function at x
%
% ----------------------------------------------------------------------
%
% Author:           Yulong Dong, dongyl@berkeley.edu
% Version:          1.0
% Last revision:    11/16/2019
%
% ----------------------------------------------------------------------

function ret = ChebyCoef2Func(x, coef, evenparity, partialcoef)
ret = zeros(length(x), 1);
y = acos(x);
if partialcoef
    if evenparity
        for k = 1:length(coef)
            ret = ret + coef(k) * cos(2*(k-1)*y);
        end
    else
        for k = 1:length(coef)
            ret = ret + coef(k) * cos((2*k-1)*y);
        end
    end
else
    if evenparity
        for k = 1:2:length(coef)
            ret = ret + coef(k) * cos((k-1)*y);
        end
    else
        for k = 2:2:length(coef)
            ret = ret + coef(k) * cos((k-1)*y);
        end
    end
end
end