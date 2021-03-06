function PKL=PeakLock(varargin)
%Determines the peak lock
%   This function determines the Peak Lock of the PIV measurement. It takes
%   a velocity component (matrix or cell) and it computes the Peak Lock
%   using a number of bins given by the user. If not given 20 bins are
%   used.
if nargin==1
    u=varargin{1};
    nbins=20;
    PKL=PeakLockCalc(u,nbins);
elseif nargin==2
    [u,nbins]=deal(varargin{:});
    PKL=PeakLockCalc(u,nbins);
elseif nargin==3
    [u,v,nbins]=deal(varargin{:});
    PeakLockCalc(u,nbins,flag);
    PKL=PeakLockCalc(v,nbins,flag);
end




