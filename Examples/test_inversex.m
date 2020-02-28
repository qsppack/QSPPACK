function test_inversex
%--------------------------------------------------------------------------
% # Test case 3: Matrix inversion
%
% We would like to approximate 1/x over [1/kappa,1] by a polynomial, such polynomial is generated
% by Remez algorithm and the approximation error is bounded by 10^{-14}
%
% parameters
%     kappa: parameters of approximation
%     parity: parity of approximation polynomial (0 -- even, 1 -- odd)
%     criteria: stop criteria, default 1e-12
%     plot_phase: whether plot phase factors
%
% note that coefficients corresponds to parameter kappa are pre-calculated,
% one can only choose kappa = 10,20,30,40,50
%
%--------------------------------------------------------------------------
% setup parameters

kappa = 50;
parity = 1;
criteria = 1e-12;
plot_phase = true;

%--------------------------------------------------------------------------
% find phase factors

opts.criteria = criteria;
if(parity==0)
    load("coef_xeven_"+int2str(kappa)+"_"+int2str(14)+".mat","coef");
else
    load("coef_xodd_"+int2str(kappa)+"_"+int2str(14)+".mat","coef");
end
if(parity==0)
    maxorder = length(coef)*2-2;
else
    maxorder = length(coef)*2-1;
end
[phi,out] = QSP_solver(coef/2,parity,opts);
parity_label = ["even" "odd"];
fprintf("- Info: \t\tQSP phase factors --- solved by L-BFGS\n")
fprintf("- Parity: \t\t%s\n- Degree: \t\t%d\n", parity_label(parity+1), maxorder);
fprintf("- Iteration times: \t%d\n", out.iter);
fprintf("- CPU time: \t%.1f s\n", out.time);

%--------------------------------------------------------------------------
% plot phase factors

if(plot_phase)
    figure(1);
    plot(1:length(phi),phi);
end

%--------------------------------------------------------------------------

end