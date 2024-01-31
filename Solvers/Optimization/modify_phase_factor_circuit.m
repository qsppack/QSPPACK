function phi_mod = modify_phase_factor_circuit(phi_full)
%--------------------------------------------------------------------------
% Construct the modified phase factor for qsvt circuit. 
%
% It preserves the symmetry when phi_full is symmetric.
%--------------------------------------------------------------------------
% Author:     Lin Lin   update 01/2024
%
%--------------------------------------------------------------------------

dd = length(phi_full)-1;
if( dd == 0 );
  phi_mod = phi_full;
  return;
end

phi_mod = zeros(size(phi_full));
phi_mod(2:end-1) = phi_full(2:end-1) + pi/2;

switch mod(dd, 4)
  case 0
    fac = pi/4;
  case 1
    fac = 0;
  case 2
    fac = -pi/4;
  case 3
    fac = pi/2;
end

phi_mod(1) = phi_full(1) + fac;
phi_mod(end) = phi_full(end) + fac;
end
