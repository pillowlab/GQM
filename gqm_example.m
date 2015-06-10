
%% 1) Add zaso to load path; set global constants
addpath('./zaso');
verbose = false;
isPoisson = false;

nT = 1e4; % number of timebins to simulate
nx = 10;  % dimensionality of stimulus in space 
          % number of time bins relevant to each response 
nd = 5;   % (the effective stimulus dimensionality at each time point will be dx = nx*nd)

dx = nx*nd;

%% 2) Form a ground truth gtModel in GQM family (canonical parameterization)
gtModel = struct();
gtModel.a = -3;
gtModel.bx = 0.1 * (2 * (1:dx)' - 8);
gtModel.Wxx = cos(0.2*pi+(1:dx)'*5/dx)';
gtModel.Dxx = -1;

if isPoisson
    gtModel.inverseLink = gqm_inverse_link_function_factory('Poisson');
else
    gtModel.inverseLink = gqm_inverse_link_function_factory('gaussian');
end

%% 3) Generate stimulus and simulate GQM response
x0 = randn(1e4, nx); % stimulus (in dimensions of (TIME x SPACE)
x = makeStimRows(x0, nd);
y = gtModel.inverseLink(gqm_evaluate_Q(gtModel, x)); % Response variable

% Simulate noisy response
if isPoisson
    y = poissrnd(y);
else
    y = y + randn(size(y)) * 1.2;
end

%% 4) Pack the stimulus into a zaso
zaso = encapsulateRaw(x, y(:));

%% 5) Estimate moments
moments = gqm_compute_moments(zaso);

%% 6) Set-up data structure (single subunit)
data(1).zaso = zaso;
data(1).moments = moments;

%% 7) Initialize subunit
if isPoisson
    model = gqm_init_model_from_moments(moments, 'poisson-gaussian');
else
    model = gqm_init_model_from_moments(moments, 'gaussian-whitened-orthant-symmetric');
end

%% 8) Predict training response with learned model
% the 'model' struct contains the fit GQM
ypred = gqm_evaluate_Q(model, x); % for Poisson noise model, ypred is the predicted rate.

%% 8) Plot Results

figure(1); clf; 
subplot(2,2,1:2)
hold on
% plot only the first 250 samples
plot(y(1:250), 'k', 'linewidth', 2) 
plot(ypred(1:250), 'r', 'linewidth', 2)
grid on
axis tight 
legend('data', 'gqm fit')
xlabel('time bin')
ylabel('response magnitude')
subplot(2,2,3)
imagesc(model.Cxx)
axis equal
axis tight
title('MEL estimate of C')

%% 9) Visualize the receptive fields by the eigendecomposition of C

[u, e] = svd(model.Cxx); 

subplot(2,2,4)
imagesc(reshape(u(:,1), nd, nx)) 
axis equal
axis tight
title('most significant receptive field')
