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
%
% Reference: Yulong Dong, Xiang  Meng, K.Birgitta Whaley and Lin Lin
%            Efficient Phase Factor Evaluation in Quantum Signal Processing
%
% Author: Yulong Dong, Xiang Meng
% Version 1.0
% Last Update 06/2020
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
    temp = phi;
    temp(1) = temp(1) - pi/4;
    temp(end) = temp(end) - pi/4;
    plot(1:length(temp),temp);
    
    x = linspace(0,1,1000);
    targ = @(x) ChebyCoef2Func(x,coef,parity,true);
    y = zeros(size(x));
    yqsp = zeros(size(x));
    for jj = 1:length(x)
        y(jj) = targ(x(jj));
        yqsp(jj) = QSPGetUnitary(phi, x(jj));
    end
    scale_fac = mean(rmmissing(yqsp./y));
    fprintf("- Linf approximation error: \t%.2e\n", norm(y*scale_fac - yqsp, inf));
end

%--------------------------------------------------------------------------
