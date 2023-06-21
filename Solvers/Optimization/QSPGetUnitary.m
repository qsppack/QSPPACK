% Get QSP unitary matrix based on given phase vector and point x \in [-1,
% 1]
%
% ----------------------------------------------------------------------
%
% Input:
%       phase, x
% Output:
%       targ    -- QSP approximation of target, real(ret(1, 1))
%
% ----------------------------------------------------------------------
% Author:           Yulong Dong      update 3/11/2020
%
%  ----------------------------------------------------------------------

function targ = QSPGetUnitary(phase, x)

    Wx = [x, 1j*sqrt(1-x^2); 1j*sqrt(1-x^2), x];
    expphi = exp(1j*phase);

    ret = [expphi(1), 0; 0, conj(expphi(1))];

    for k = 2:numel(expphi)
        temp = [expphi(k), 0; 0, conj(expphi(k))];
        ret = ret * Wx * temp;
    end

    targ = real(ret(1,1));

end

