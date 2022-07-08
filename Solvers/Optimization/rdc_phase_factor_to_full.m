function [phi_full] = rdc_phase_factor_to_full(phi_cm, parity,targetPre)
%--------------------------------------------------------------------------
% construct the full phase factors from reduced phase factors             
%--------------------------------------------------------------------------

phi_right = phi_cm;
if (targetPre==true)
    phi_right(end) = phi_right(end) + pi/4;
end
dd = 2 * length(phi_right);
if parity == 0
    dd = dd - 1;
end
phi_full = zeros(dd, 1);
phi_full((dd+1-length(phi_right)):end) = phi_right;
phi_full(1:length(phi_right)) = phi_full(1:length(phi_right)) + phi_right(end:-1:1);

end