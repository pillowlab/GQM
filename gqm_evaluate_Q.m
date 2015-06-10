function [q y dy ddy] = gqm_evaluate_Q(model, x)
% [q y dy ddy] = gqm_evaluate_Q(model, x);
% Evaluate the GQM model in canonical parameters.
%
% Input
%   model: model structure with inverseLink and canonical parameterization
%   x: [N x dx] stimulus matrix
%
% Warning: This is an internal function that should be called through zaso only.
%
% Compute quadratic function of a given stimulus vector x. 
%   q(x) = a + b'*x + x'*C*x
%
% and the output
%   y(x) = f(q(x));
%
% Expects model struct containing parameters a, bx, Wxx, and Dxx, 
% where Wxx and Dxx are the eigenvectors and vector of eigenvalues 
% of C, respectively. 

assert(isa(model, 'struct'));

q = model.a;

if isfield(model, 'bx') && ~isempty(model.bx)
    q = q + model.bx'*x';
end

if isfield(model, 'Wxx') && ~isempty(model.Wxx) % canonical parameterization
    q = q + model.Dxx'*((model.Wxx*x').^2) ;
elseif isfield(model, 'Cxx') && ~isempty(model.Cxx)
    q = model.a + model.bx'*x' + sum(bsxfun(@times, (model.Cxx*x'),x'));
else
    error('Invalid model structure.')
end

if nargout > 1
    [y dy ddy] = model.inverseLink(q);
end
