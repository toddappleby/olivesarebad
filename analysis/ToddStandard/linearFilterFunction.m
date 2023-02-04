function f = linearFilterFunction(theta, t)
% f = linearFilterFunction(theta, t)

numBasis = length(theta);

if numBasis == 5
    % Define the temporal filter function.
    tfun = @(p,t)(p(1) .* (((t./p(2)).^4)./(1+((t./p(2)).^4))) .* exp(-((t./p(3)))) .* cos(((2.*pi.*t)./p(4))+(2*pi*p(5)/360)));
    f = tfun(theta, t);
elseif numBasis == 6
    tfun = @(p,t)(p(1) .* (((t./p(2)).^p(6))./(1+((t./p(2)).^p(6)))) .* exp(-((t./p(3)))) .* cos(((2.*pi.*t)./p(4))+(2*pi*p(5)/360)));
    f = tfun(theta, t);
else
    t2 = 2*t - t.^2;
    m1 = t2(:) * ones(1,15);
    m2 = ones(length(t),1) * (1:15);
    A = sin(m1 .* m2 * pi);
    basis = orth(A);
    basis = basis(:,numBasis:end);
    f = basis * theta'; % The filter
    f = f(:)';
end