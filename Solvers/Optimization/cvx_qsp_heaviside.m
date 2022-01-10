function phi_qc = cvx_qsp_unwrap(deg, mu, sigma_min, sigma_mu_m, ...
  sigma_mu_p, sigma_max, npts, epsil, fscale, criteria)
%% Return the phase factors for a step function
%
% The step function takes the form
%     f(x) = 0,    x <= mu
%            1,    x >  mu
%
% The purpose of this function is to interface with python.
%
% Lin Lin
% Last revision: 1/6/2022

% This is due to that the default class of numerical integers in
% matlab is "double"

deg = double(deg);
npts = double(npts);

opts.npts = npts; 
opts.epsil = epsil;
opts.fscale = fscale;
opts.criteria = criteria;
opts.intervals= [sigma_min, sigma_mu_m, sigma_mu_p, sigma_max];
opts.isplot = false;
% opts.fname = 'func';
opts.objnorm = inf;

func = @(x) (x > mu);
phi_qc = cvx_qsp(func, double(deg), opts);
