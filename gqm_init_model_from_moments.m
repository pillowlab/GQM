function model = gqm_init_model_from_moments(moments, likelihood_method)
% Fit the GQM model without history terms by moment based expected ML solutions.
% Under the assumption that the stimulus is Gaussian, and the canonical link
% function is used, a fast moment based estimation of model parameters can be
% obtained via expected-likelihood trick.
%
% Moments up to fourth order may be required to be finite.
%
% TODO take prior distribution of input
% TODO take prior as regularizer
% TODO regularize Poisson-gaussian case
%
% If whitened signal is used for the estimation, the original solution can be
% obtained by multiplying the whitening matrix U to get U*bx and U*Cxx*U'
%
% Input
%   moments: structure obtained from gqm_moments_timeseries.m
%	    If the exact stimulus distribution is known replace corresponding
%	    entries in the moments structure.
%   likelihood_method: string option
%	'gaussian-gaussian': Gaussian likelihood with Gaussian stimuli
%	'gaussian-whitened-sym-id': Gaussian likelihood with whitened stimulus
%	    (EXX = I) distribution with identically distributed symmetric 
%	    marginals (mean zero P(x) = P(-x)).
%	'gaussian-whitened-orthant-symmetric': Gaussian likelihood with
%	    whitened stimulus that is orthant symmetric, i.e.,
%	    P(x_1, x_2, ...) = P(rho_1 x_1, rho_2 x_2, ...) for any
%	    combination of rho_i = (+1 or -1).
%	'poisson-gaussian': Poisson likelihood with gaussian stimulus.
%
% Output is a single structure
%   model.a: (1) bias
%   model.bx: (p*d x 1) linear filter
%   model.Cxx: (p*d x p*d) quadratic filter
%   model.method: likelihood_method string
%   model.Id: version information string
%
% $Id$

switch(lower(likelihood_method))
case 'gaussian-gaussian'
    % same as (Koh and Powers, 1985)
    bx  = moments.EXX \ moments.EYX;
    Cxx = (moments.EXX \ moments.EYXX) / moments.EXX / 2;
    a   = -trace(Cxx * moments.EXX) + moments.EY;
    model.likelihood = 'Gaussian';
case 'gaussian-whitened-sym-id'
    % based on (Wiens, 1992)
    bx = moments.EYX;
    mu4 = sum(diag(moments.EXXXX)); % diagonal terms should be constant
    mu22 = sum(moments.EXXXX(:)) - mu4; % off-diagonal terms should be constant
    mu4 = mu4 / moments.dim; mu22 = mu22 / (moments.dim^2 - moments.dim);
    Cxx_diag = diag((diag(moments.EYXX) - moments.EY) / (mu4 - mu22));
    Cxx_offdiag = (moments.EYXX/2) ./ mu22;
    Cxx_offdiag = Cxx_offdiag - diag(diag(Cxx_offdiag));
    Cxx =  Cxx_diag + Cxx_offdiag;
    a = -trace(Cxx) + moments.EY;
    model.likelihood = 'Gaussian';
case 'gaussian-whitened-orthant-symmetric'
    % Evan and Memming's derivation for non-identical marginal distributions
    bx = moments.EYX;
    Cxx_diag = diag((diag(moments.EYXX) - moments.EY)' / (moments.EXXXX - 1));
    Cxx_offdiag = (moments.EYXX/2) ./ moments.EXXXX;
    Cxx_offdiag = Cxx_offdiag - diag(diag(Cxx_offdiag));
    Cxx =  Cxx_diag + Cxx_offdiag;
    a = -trace(Cxx) + moments.EY;
    model.likelihood = 'Gaussian';
case 'poisson-gaussian'
    % very similar to (Park & Pillow 2011)
    % involves two matrix inversions
    assert(moments.EY >= 0, 'non-negative mean spike count is a must')
    LambdaPrime = moments.EYXX + moments.EYX * moments.EYX' * 2 / moments.EY;
    %LambdaPrime = LambdaPrime + 1e-5 * eye(moments.dim); % regularizer
    %invLambdaPrime = pinv(LambdaPrime);
    bx = LambdaPrime \ moments.EYX;
    %bx = pinv(LambdaPrime + 1e-5 * eye(moments.dim)) * moments.EYX;
    %invEXX = pinv(moments.EXX);
    invEXX = inv(moments.EXX);
    Cxx = (invEXX - moments.EY * pinv(LambdaPrime)) / 2;
    a = moments.EY * sqrt(det(eye(size(moments.EXX)) - 2 * Cxx * moments.EXX)) + 0.5 * bx' * ((invEXX - 2 * Cxx) \ bx);
    if ~isreal(a) || isnan(a)
	fprintf('Cxx estimate is likely unstable...trying to fix it [%s]\n', mfilename);
	Cxx = (invEXX - moments.EY * pinv(1e-5 * eye(moments.dim) + LambdaPrime)) / 2;
	a = moments.EY * sqrt(det(eye(size(moments.EXX)) - 2 * Cxx * moments.EXX)) + 0.5 * bx' * ((invEXX - 2 * Cxx) \ bx);
	if ~isreal(a) || isnan(a)
	    fprintf('Cxx estimate is no no no ...trying to fix it [%s]\n', mfilename);
	    a = moments.EY * 0.9;
	end
    end
    model.likelihood = 'Poisson';
otherwise
    error('Unknown likelihood-method [%s]', likelihood_method);
end

model.a = a;
model.bx = bx;
model.Cxx = Cxx;
model.fit_ML_method = likelihood_method;
model.Id = '$Id$';
