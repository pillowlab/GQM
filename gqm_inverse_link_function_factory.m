function f = gqm_inverse_link_function_factory(name)
% f = gqm_inverse_link_function_factory(name);
% returns the point nonlinearity (aka mean function or inverse link function).
%
% Input
%   name: could be the name of the nonlinear function or a exponential family
%	distribution. In the latter case, it returns the CANONICAL inverse link.
%
% Output
%   f: @(x)->[f df ddf]
%
% Note: more functions could be found in "ncclabcode/nlfuns"

switch lower(name)
    case {'exponential', 'poisson'}
	f = @expfun;
    case 'logexp1'
	f = @logexp1;
    case {'identity', 'gaussian'}
        f = @identity;
    otherwise
        error('Unknown link function [%s]', name);
end

end % gqm_inverse_link_function_factory

function [f df ddf] = identity(x)
    f   = x;
    df  = ones(size(x));
    ddf = zeros(size(x));
end

function [f,df,ddf] = expfun(x)
    f = exp(x);
    df = f;
    ddf = df;
end

function [f,df,ddf] = logexp1(x)
    % [f,df,ddf] = logexp1(x);
    %
    % Computes the function:
    %    f(x) = log(1+exp(x))
    %
    % and returns first and second derivatives

    f = log(1+exp(x));

    if nargout > 1
	df = exp(x)./(1+exp(x));
    end

    if nargout > 2
	ddf = exp(x)./(1+exp(x)).^2;
    end

    % Check for small values
    if any(x(:)<-20)
	iix = (x(:)<-20);
	f(iix) = exp(x(iix));
	df(iix) = f(iix);
	ddf(iix) = f(iix);
    end

    % Check for large values
    if any(x(:)>500)
	iix = (x(:)>500);
	f(iix) = x(iix);
	df(iix) = 1;
	ddf(iix) = 0;
    end
end
