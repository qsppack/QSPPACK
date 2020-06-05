%--------------------------------------------------------------------------
% Test case 1: Hamiltonian simulation
%
% Here we want to approxiamte e^{-i\tau x} by Jacobi-Anger expansion:
% 
% e^{-i\tau x} = J_0(\tau)+2\sum_{k even} (-1)^{k/2}J_{k}(\tau)T_k(x)+2i\sum_{k odd} (-1)^{(k-1)/2}J_{k}(\tau) T_k(x)
%
% We truncate the series up to N = 1.4\tau+log(10^{14}), which gives an polynomial approximation of e^{-i\tau x} with
% accuracy 10^{-14}. Besides, we deal with real and imaginary part of the truncated series seperatly and divide them
% by a constant factor 2 to enhance stability.
%
% parameters
%     tau: the duration \tau in Hamiltonian simulation
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

tau = 1000;
criteria = 1e-12;
plot_phase = true;

%--------------------------------------------------------------------------
% find phase factors

opts.criteria = criteria;
maxorder = ceil(1.4*tau+log(1e14));
if(mod(maxorder,2)==1)
    maxorder = maxorder-1;
end

%--------------------------------------------------------------------------
% even part

coef = zeros(maxorder/2+1,1);
for i=1:length(coef)
    coef(i) = (-1)^(i-1)*besselj(2*(i-1),tau);
end
coef(1) = coef(1)/2;
[phi1,out1] = QSP_solver(coef,0,opts);

fprintf("- Info: \t\tQSP phase factors --- solved by L-BFGS\n")
fprintf("- Parity: \t\t%s\n- Degree: \t\t%d\n", "even", maxorder);
fprintf("- Iteration times: \t%d\n", out1.iter);
fprintf("- CPU time: \t%.1f s\n", out1.time);

%--------------------------------------------------------------------------
% odd part

coef = zeros(maxorder/2+1,1);
for i=1:length(coef)
    coef(i) = (-1)^(i-1)*besselj(2*i-1,tau);
end
[phi2,out2] = QSP_solver(coef,1,opts);

%--------------------------------------------------------------------------
% output

fprintf("- Info: \t\tQSP phase factors --- solved by L-BFGS\n")
fprintf("- Parity: \t\t%s\n- Degree: \t\t%d\n", "odd", maxorder+1);
fprintf("- Iteration times: \t%d\n", out2.iter);
fprintf("- CPU time: \t%.1f s\n", out2.time);

%--------------------------------------------------------------------------
% plot phase factors

if(plot_phase)
    figure(1);
    temp = phi1;
    temp(1) = temp(1) - pi/4;
    temp(end) = temp(end) - pi/4;
    plot(1:length(temp),temp);
    
    x = linspace(0,1,1000);
    y = cos(tau*x);
    yqsp = zeros(size(x));
    for jj = 1:length(x)
        yqsp(jj) = QSPGetUnitary(phi1, x(jj));
    end
    scale_fac = mean(rmmissing(yqsp./y));
    fprintf("- Linf approximation error (even): \t%.2e\n", norm(y*scale_fac - yqsp, inf));
    
    figure(2);
    temp = phi2;
    temp(1) = temp(1) - pi/4;
    temp(end) = temp(end) - pi/4;
    plot(1:length(temp),temp);
    
    x = linspace(0,1,1000);
    y = sin(tau*x);
    yqsp = zeros(size(x));
    for jj = 1:length(x)
        yqsp(jj) = QSPGetUnitary(phi2, x(jj));
    end
    scale_fac = mean(rmmissing(yqsp./y));
    fprintf("- Linf approximation error (odd): \t%.2e\n", norm(y*scale_fac - yqsp, inf));
end

%--------------------------------------------------------------------------
