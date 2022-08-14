%--------------------------------------------------------------------------
% Test case : Chebyshev polynomial
%
% Generate the phase factors for the Chebyshev polynomial that can be
% readily used in a quantum circuit.
%
% parameters
%     d: degree of the Chebyshev polynomial
%
%--------------------------------------------------------------------------
%
% Reference: Yulong Dong, Xiang  Meng, K.Birgitta Whaley and Lin Lin
%            Efficient Phase Factor Evaluation in Quantum Signal Processing
%
%--------------------------------------------------------------------------
% setup parameters

d = 10;
parity = mod(d,2);
if parity == 0 
  coef = zeros(round(d/2+1), 1);
  coef(end) = 1.0;
else
  coef = zeros(round((d+1)/2), 1);
  coef(end) = 1.0;
end
criteria = 1e-12;

%--------------------------------------------------------------------------
% find phase factors

opts.criteria = criteria;
[phi,out] = QSP_solver(coef,parity,opts);
parity_label = ['even', 'odd'];
fprintf('- Info: \t\tQSP phase factors --- solved by L-BFGS\n')
fprintf('- Parity: \t\t%s\n- Degree: \t\t%d\n', parity_label(parity+1), d);
fprintf('- Iteration times: \t%d\n', out.iter);
fprintf('- CPU time: \t%.1f s\n', out.time);
disp(phi)

fprintf('Generate the phase factors that can be directly used in the quantum circuit\n');
phi_qc = zeros(size(phi));
phi_qc(1) = phi(1)+pi/4;
phi_qc(2:end-1) = phi(2:end-1)+pi/2;
phi_qc(end) = phi(end)+pi/4;
disp(phi_qc)
