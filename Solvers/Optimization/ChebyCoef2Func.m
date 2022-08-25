% Evaluate function based on Chebyshev expansion coefficients
%
% ----------------------------------------------------------------------
%
% Input:
%       x, coef
%       parity  -- 0 for even, 1 for odd
%  partialcoef  -- true: only include even/odd coefficiennts
% Output:
%       ret     -- function at x
%
% ----------------------------------------------------------------------
% Author:           Yulong Dong  update 06/2020
%
% ----------------------------------------------------------------------

function ret = ChebyCoef2Func(x, coef, parity, partialcoef)
ret = zeros(length(x), 1);
y = acos(x);
if partialcoef
    if parity == 0
        for k = 1:length(coef)
            ret = ret + coef(k) * cos(2*(k-1)*y);
        end
    else
        for k = 1:length(coef)
            ret = ret + coef(k) * cos((2*k-1)*y);
        end
    end
else
    if parity == 0
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
