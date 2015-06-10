function moments = gqm_comptue_moments(zaso)
% moments = gqm_comptue_moments(zaso);
% Estimate moments from a zaso.
%
% Input:
%   zaso: structure that spits out data when poked
%      zaso.N: (1) number of samples that it will give us
%      zaso.X: @(n) -> (1 x p) returns input vector
%      zaso.Y: @(n) -> (1) returns output
%
% Output is a single structure
%   moments.EY: (1) mean of output E[Y(t)]
%   moments.EX: (p*d x 1) mean of stimulus E[X(t)]
%   moments.EYX: (p*d x 1) cross-correlation E[Y(t) X(t)]
%   moments.EXX: (p*d x p*d) auto-correlation E[X(t) X(t)']
%   moments.EYXX: (p*d x p*d) cross-bicorrelation E[Y(t) X(t) X(t)']
%   moments.EXXXX: (p*d x p*d) 4th order moments E[X_i(t)^2 X_j(t)^2]
%   moments.nAvg: (1) number of samples used to estimate the moments
%   moments.Id: version information string
%
% $Id$

rsum = zasoFarray(zaso, {@(x,y) sum(x)', @(x,y) sum(y), @(x,y) x' * y, @(x,y) x' * x, @(x,y) bsxfun(@times, x, y)' * x, @(x,y) (x.^2)' * (x.^2)}, {});

moments.EX = rsum{1} / zaso.N;
moments.EY = rsum{2} / zaso.N;
moments.EYX = rsum{3} / zaso.N;
moments.EXX = rsum{4} / zaso.N;
moments.EYXX = rsum{5} / zaso.N;
moments.EXXXX = rsum{6} / zaso.N;
moments.dim = size(moments.EX, 1);
moments.Id = '$Id$';
