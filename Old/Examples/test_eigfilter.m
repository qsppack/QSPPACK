%--------------------------------------------------------------------------
% Test case 2: Eigenstate filter
%
% Here we want to generate factors for the eigenstate filter function:
% 
% f_n(x,\delta)=\frac{T_n(-1+2\frac{x^2-\delta^2}{1-\delta^2})}{T_n(-1+2\frac{-\delta^2}{1-\delta^2})}.
%
% We divide f_n by a constant factor \sqrt{2} to enhance stability.
%
% Reference: Lin Lin and Yu Tong
%            Solving quantum linear system problem with near-optimal complexity
%
% parameters
%     n,delta: parameters of f_n
%     criteria: stop criteria, default 1e-12
%     plot_phase: whether plot phase factors
%
% note that coefficients corresponds to pair (n,delta) are pre-calculated,
% one can only choose delta=0.1,0.05,0.01,0.005 and 
% n*delta = 3,5,10,15,20,25
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

n = 1000;
delta = 0.01;
criteria = 1e-12;
plot_phase = true;

%--------------------------------------------------------------------------
% find phase factors

opts.criteria = criteria;
maxorder = 2*n;
load("coef_"+int2str(n)+"_"+num2str(delta)+".mat","coef");
[phi,out] = QSP_solver(coef,0,opts);

fprintf("- Info: \t\tQSP phase factors --- solved by L-BFGS\n")
fprintf("- Parity: \t\t%s\n- Degree: \t\t%d\n", "even", maxorder);
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
    temp_fac = cos(n*acos(-1-2*delta^2/(1-delta^2)));
    y = cos(n*acos(-1+2*(x.^2-delta^2)/(1-delta^2)));
    y = real(y / temp_fac);
    yqsp = zeros(size(x));
    for jj = 1:length(x)
        yqsp(jj) = QSPGetUnitary(phi, x(jj));
    end
    scale_fac = mean(rmmissing(yqsp./y));
    fprintf("- Linf approximation error: \t%.2e\n", norm(y*scale_fac - yqsp, inf));
end

%--------------------------------------------------------------------------
