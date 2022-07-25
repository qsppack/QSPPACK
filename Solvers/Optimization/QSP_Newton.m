function [phi, err, iter, runtime] = QSP_Newton(coef,parity,opts)
%--------------------------------------------------------------------------
% Newton solver for finding phase factors such that the real part of the 
% (1,1) element of the QSP unitary matrix gives desire Chebyshev expansion.

% Input:
%         coef --- Chebyshev coefficients
%       parity --- Parity of phi (0 -- even, 1 -- odd)
%          opts --- Options structure with fields
%                   maxiter: maximal iteration number
%                   criteria: stop criteria 
%                   targetPre: Pre to be target function 
%                   useReal: use real matrix mulplitcation to get QSP entry
%
% Output:
%     phi_full --- full set of phase factors
%                  
%--------------------------------------------------------------------------
% setup options for Newton solver
if ~isfield(opts,'maxiter');              opts.maxiter = 1e5; end
if ~isfield(opts,'criteria');             opts.criteria = 1e-12; end
if ~isfield(opts,'targetPre');            opts.targetPre = true;    end
if ~isfield(opts,'useReal');              opts.useReal = true; end
if ~isfield(opts,'print');                opts.print = 1; end
if ~isfield(opts,'itprint');              opts.itprint = 1; end

tic

%--------------------------------------------------------------------------
% copy value to parameters

maxiter = opts.maxiter;
crit = opts.criteria;
pri = opts.print;
itprint = opts.itprint;   

%--------------------------------------------------------------------------
%  setup print format

stra1 = ['%4s','%13s','\n'];
str_head = sprintf(stra1,'iter','err');
str_num = '%4d  %+5.4e \n';


%%--------------------------------------------------------------------------
% initial preparation
if (opts.targetPre == true)   
    coef = - coef; % inverse is necessary
end
phi = coef/2;
err = inf;
iter = 0;

%--------------------------------------------------------------------------
% solve by Newton method

while err > crit && iter < maxiter
    [Fval,DFval] = F_Jacobian(phi, parity,opts);
    res = Fval - coef;
    phi = phi - DFval\res;
    err = norm(res, 1);
    iter = iter + 1;
    if(pri&&mod(iter,itprint)==0)
        if(iter==1||mod(iter-itprint,itprint*10)==0)
            fprintf("%s",str_head);
        end
        fprintf(str_num,iter,err);
    end
    if iter>=maxiter; fprintf("Max iteration reached.\n"); break; end
    if err<crit; fprintf("Stop criteria satisfied.\n"); break; end
end

runtime = toc;
end





