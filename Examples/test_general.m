function test_general
%--------------------------------------------------------------------------
% Approximating a general polynomial expanded in Chebyshev basis. You may 
% need to generate the target polynomial first, e.g., by the Remez method
% provided in Remez.ipynb.
%
% parameters
%     data_name: name of the data file that contains
%                    coef: the coefficients of polynomial
%                    parity: the parity of the polynomial
%     criteria: stop criteria, default 1e-12
%     plot_phase: whether plot phase factors
%
%--------------------------------------------------------------------------
% setup parameters

data_name = "";
criteria = 1e-12;
plot_phase = true;

%--------------------------------------------------------------------------
% find phase factors

load(data_name,"coef","parity");
opts.criteria = criteria;
if(parity==0)
    maxorder = length(coef)*2-2;
else
    maxorder = length(coef)*2-1;
end
[phi,out] = QSP_solver(coef,parity,opts);
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