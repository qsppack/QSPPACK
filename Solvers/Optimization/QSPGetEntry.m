% Get QSP unitary matrix based on given phase vector and a list of points x
% \in [-1,1]
%
% ----------------------------------------------------------------------
%
% Input:
%       xlist, phase: full phase factors
%       opts: contains information of phase factors
%             targetPre: get Pre if true, otherwise, get Pim
%             parity, typePhi
% Output:
%       ret    -- QSP approximation of target
%
% ----------------------------------------------------------------------
%
% Author:     Jiasu Wang   update 07/2020
% 
%--------------------------------------------------------------------------

function ret= QSPGetEntry(xlist, phase, opts)

typePhi = opts.typePhi;
targetPre = opts.targetPre;
parity = opts.parity;

d = length(xlist);
ret = zeros(d,1);

if strcmp(typePhi,'reduced')
    dd = 2 * length(phase)-1+parity;
    phi = zeros(dd, 1);
    phi((dd+1-length(phase)):end) = phase;
    phi(1:length(phase)) = phi(1:length(phase)) + phase(end:-1:1);
else
    phi = phase;
end

if targetPre == false
    phi(1) = phi(1)-pi/4;
    phi(end) = phi(end)-pi/4;
end
for i = 1:d
    x = xlist(i);
    ret(i) = QSPGetUnitary(phi, x);
end
end
