% Get QSP unitary matrix based on given modified phase vector and point
% x \in [-1, 1]
%
% The modified phase factors are used in the qsvt circuit.
%
% ----------------------------------------------------------------------
%
% Input:
%       phi_mod: modified phase factor of length d+1. This is converted
%         using the modify_phase_factor_circuit routine that preserves
%         the symmetry when possible.
%       x: real number in [-1,1].
% Output:
%       targ    -- (1,1) entry of the QSP unitary
%
% ----------------------------------------------------------------------
% Author:           Lin Lin          update 01/30/2024
%  ----------------------------------------------------------------------

function targ = QSPGetUnitary_modified_phase(phi_mod, x)

    % form of Ux from the Hermitian block encoding (or CS decomposition)
    Ux = [x, sqrt(1-x^2); sqrt(1-x^2), -x];
    
    expphi = exp(1j*phi_mod);

    ret = [expphi(1), 0; 0, conj(expphi(1))];

    for k = 2:numel(expphi)
        temp = [expphi(k), 0; 0, conj(expphi(k))];
        ret = ret * Ux * temp;
    end

    targ = ret(1,1);

end

