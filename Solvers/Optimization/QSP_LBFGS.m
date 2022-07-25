function [phi,obj_value,iter] = QSP_LBFGS(obj,grad,delta,phi,opts)
%--------------------------------------------------------------------------
% Solving phase factors optimization via L-BFGS
%
% Input:
%        obj --- Objective function L(phi) (Should also be given in grad)
%       grad --- Gradient of obj function
%      delta --- Samples
%        phi --- Initial value
%       opts --- Options structure with fields
%              maxiter: max iteration 
%              gamma: linesearch retraction rate
%              accrate: linesearch accept ratio
%              minstep: minimal stepsize
%              criteria: stop criteria for obj value on Chevyshev points
%              lmem: L-BFGS memory size
%              print: whether to output
%              itprint: print frequency
%              parity: parity of polynomial (0 -- even, 1 -- odd)
%              target: target polynomial
%
% Output:
%         phi --- Solution of phase factors optimization
%   obj_value --- Objection value at optimal point L(phi^*)
%        iter --- iteration number
%
%--------------------------------------------------------------------------
%
% Reference: Yulong Dong, Xiang  Meng, K.Birgitta Whaley and Lin Lin
%            Efficient Phase Factor Evaluation in Quantum Signal Processing
%
% Version 2.0
% Last Update 07/2022
%
%--------------------------------------------------------------------------
% options for L-BFGS solver

if ~isfield(opts,'maxiter');              opts.maxiter = 5e4; end
if ~isfield(opts,'gamma');                opts.gamma = 0.5; end
if ~isfield(opts,'accrate');              opts.accrate = 1e-3; end
if ~isfield(opts,'minstep');              opts.minstep = 1e-5; end
if ~isfield(opts,'criteria');             opts.criteria = 1e-12; end
if ~isfield(opts,'lmem');                 opts.lmem = 200; end
if ~isfield(opts,'print');                opts.print = 1; end
if ~isfield(opts,'itprint');              opts.itprint = 1; end

%--------------------------------------------------------------------------
% copy value to parameters

maxiter = opts.maxiter;   gamma = opts.gamma;      accrate = opts.accrate;
lmem = opts.lmem;         minstep = opts.minstep;  pri = opts.print;
itprint = opts.itprint;   crit = opts.criteria;    

%--------------------------------------------------------------------------
%  setup print format

stra1 = ['%4s','%13s','%10s','%10s','\n'];
str_head = sprintf(stra1, ...
    'iter','obj','stepsize','des_ratio');
str_num = '%4d  %+5.4e %+3.2e %+3.2e\n';

%--------------------------------------------------------------------------
% initial computation

iter = 0;
d = length(phi);
mem_size = 0;
mem_now = 0;
mem_grad = zeros(lmem,d);
mem_obj = zeros(lmem,d);
mem_dot = zeros(lmem);
[grad_s, obj_s] = grad(phi,delta,opts);
obj_value = mean(obj_s);
GRAD = mean(grad_s)';

%--------------------------------------------------------------------------
% start L-BFGS algorithm

if(pri)
    fprintf('L-BFGS solver started \n');
end
while(true)
    iter = iter+1;
    theta_d = GRAD;
    alpha = zeros(mem_size,1);
    for i=1:mem_size
        subsc = mod(mem_now-i,lmem)+1;
        alpha(i) = mem_dot(subsc)*(mem_obj(subsc,:)*theta_d);
        theta_d = theta_d-alpha(i)*mem_grad(subsc,:)';
    end
    theta_d = 0.5*theta_d;
    if(opts.parity==0)
        theta_d(1) = theta_d(1)*2;
    end
    for i=1:mem_size
        subsc = mod(mem_now-(mem_size-i)-1,lmem)+1;
        beta = mem_dot(subsc)*(mem_grad(subsc,:)*theta_d);
        theta_d = theta_d+(alpha(mem_size-i+1)-beta)*mem_obj(subsc,:)';
    end
    step = 1;
    exp_des = GRAD'*theta_d;
    while(true)
        theta_new = phi-step*theta_d;
        obj_snew = obj(theta_new,delta,opts);
        obj_valuenew = mean(obj_snew);
        ad = (obj_value-obj_valuenew);
        if(ad>exp_des*accrate*step||step<minstep)
            break;
        end
        step = step*gamma;
    end
    phi = theta_new;
    obj_value = obj_valuenew;
    obj_max = max(obj_snew);
    [grad_s,~]= grad(phi,delta,opts);
    GRAD_new = mean(grad_s)';
    mem_size = min(lmem,mem_size+1);
    mem_now = mod(mem_now,lmem)+1;
    mem_grad(mem_now,:) = GRAD_new-GRAD;
    mem_obj(mem_now,:) = -step*theta_d;
    mem_dot(mem_now) = 1/(mem_grad(mem_now,:)*mem_obj(mem_now,:)');
    GRAD = GRAD_new;
    if(pri&&mod(iter,itprint)==0)
        if(iter==1||mod(iter-itprint,itprint*10)==0)
            fprintf("%s",str_head);
        end
        fprintf(str_num,iter,obj_max,step,ad/(exp_des*step));
    end
    if iter>=maxiter; fprintf("Max iteration reached.\n"); break; end
    if obj_max<crit^2; fprintf("Stop criteria satisfied.\n"); break; end
end

end
