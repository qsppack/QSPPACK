function plot_phase_factor_decay(Phi,coef,opts)
%--------------------------------------------------------------------------
% Plot the phase factor to show decay property     
%--------------------------------------------------------------------------
%
% Input:
%       phi, coef, 
%       opts.typePhi---Phi is either 'full' or 'reduced'
%       opts.targetPre---wheather Phi is added with pi/4 at the ends
% Output:
%       h --- Chebyshev coefficients up to maxorder
%
% ----------------------------------------------------------------------
if strcmp(opts.typePhi,'full')
    d = ceil((length(Phi)+1)/2);
    phi = Phi(d:end);
else
    phi = Phi;
end
if opts.targetPre == true
    phi(end)= phi(end)-pi/4;
end

figure()
plot(1:length(phi), abs(phi),'b-', 1:length(coef), abs(coef), 'r:','LineWidth',1.5);
xlabel("index",'FontSize', 15)
ylabel("magnitude",'FontSize', 15)
xlim([1 length(phi)])
legend("$\Phi$","$c$","Interpreter","latex","FontSize",15);

end