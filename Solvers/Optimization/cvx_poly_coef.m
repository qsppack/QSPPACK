function [coef_full] = cvx_poly_coef(func, deg, opts)
%% Find a polynomial approximation for a given function.
%--------------------------------------------------------------------------
%
% Input:
%       func, deg, opts
% Output:
%       c --- Chebyshev coefficients of best approximation polynomial
%
% ----------------------------------------------------------------------
% Author:           Lin Lin      update 06/2022
%                   Jiasu Wang   update 08/2022
%
% ----------------------------------------------------------------------


if( ~isfield(opts, 'npts') )
  opts.npts = 200;
end
if( ~isfield(opts, 'epsil') )
  opts.epsil = 0.01;
end
if( ~isfield(opts, 'fscale') )
  opts.fscale = 1-opts.epsil;
end
if( ~isfield(opts, 'intervals') )
  opts.intervals = [0,1];
end
if( ~isfield(opts, 'isplot') )
  opts.isplot = false;
end
if( ~isfield(opts, 'objnorm') )
  opts.objnorm = inf;
end

% check variables and assign local variables
assert(mod(length(opts.intervals), 2) == 0);
parity = mod(deg, 2);
epsil = opts.epsil;
npts = opts.npts;

[xpts, ~] = chebpts(2*npts);
xpts=union(xpts,opts.intervals');
xpts = xpts(xpts>=0);
npts = length(xpts);

n_interval = length(opts.intervals) / 2;
ind_union = [];
for i = 1 : n_interval
  ind_set{i} = find((xpts >= opts.intervals(2*i-1)) &...
    (xpts<=opts.intervals(2*i)));

  ind_union = union(ind_union, ind_set{i});
end


% evaluate the target function
fx = zeros(npts,1);
fx(ind_union) = opts.fscale * func(xpts(ind_union));


% prepare the Chebyshev polynomials
if( parity == 0 )
  n_coef = deg / 2 + 1;
else
  n_coef = (deg + 1) / 2;
end

Ax = zeros(npts, n_coef);
for k = 1 : n_coef
  if( parity == 0 )
    Tcheb = chebpoly(2*(k-1));
  else
    Tcheb = chebpoly(2*k-1);
  end
  Ax(:,k) = Tcheb(xpts);
end


%% Use CVX to optimize the Chebyshev coefficients
cvx_begin quiet
variable coef(n_coef)
variable y(npts)
minimize ( norm(y(ind_union)-fx(ind_union), opts.objnorm)) 

y==Ax*coef
y>=-(1-epsil)
y<=(1-epsil)
cvx_end

err_inf = norm(y(ind_union)-fx(ind_union), opts.objnorm);
fprintf('norm error = %g\n', err_inf);


% use chebfun to make sure the maximum is less than 1
coef_full = zeros(deg+1,1);
if( parity == 0 )
  coef_full(1:2:end) = coef;
else
  coef_full(2:2:end) = coef;
end
sol_cheb = chebfun(coef_full, 'coeffs');
max_sol = max(abs(sol_cheb));
fprintf('max of solution = %g\n', max_sol);
if( max_sol > 1.0-1e-10 )
  error('Solution is not bounded by 1. Increase npts');
end

if( opts.isplot )
  figure(1)
  clf
  hold on
  plot(xpts, y, 'ro', 'LineWidth',1.5);
  for i = 1 : n_interval
    plot(xpts(ind_set{i}), y(ind_set{i}),'b-','LineWidth',2)
  end
  hold off
  xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 15);
  ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 15);
  legend({'polynomial','target'},'FontSize', 15);


  figure(2)
  clf
  hold on
  for i = 1 : n_interval
    plot(xpts(ind_set{i}), abs(y(ind_set{i})-fx(ind_set{i})), 'k-',...
        'LineWidth',1.5)
  end
  hold off
  xlabel('$x$', 'Interpreter', 'latex','FontSize', 15);
  ylabel('$|f_\mathrm{poly}(x)-f(x)|$', 'Interpreter', 'latex','FontSize', 15);
end


end
