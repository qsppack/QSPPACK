function z = convfft(x,y,varargin)
%% Compute the convolution of two column vectors using FFT.
%--------------------------------------------------------------------------
%
% Input:
%       x, y --- two column vectors
%       varargin --- controls whether the conv is truncated
% Output:
%       z --- The convolution of x and y
%
% ----------------------------------------------------------------------
% Author:           Hongkang Ni   update 06/2025
%
% ----------------------------------------------------------------------
lenx = length(x);
leny = length(y);
truncate = false;
if ~isempty(varargin)
    truncate = varargin{1};
end
if min(lenx, leny) < 3000
    z = conv(x, y);
elseif (~truncate)
    n = 2^nextpow2(lenx+leny);
    x = [x; zeros(n-lenx,1)];
    y = [y; zeros(n-leny,1)];
    z = fft(ifft(x).* ifft(y))*n;
    z = z(1:lenx+leny-1);
else
    n = 2^nextpow2(max(lenx,leny));
    x = [x; zeros(n-lenx,1)];
    y = [y; zeros(n-leny,1)];
    z = fft(ifft(x).* ifft(y))*n;
end

if truncate
    z = z(lenx:leny);
end
    
    
    