function Y = multiHZNL(params, X)

Y(:,1) = params(1)*normcdf(params(2)*X(:,1)+params(3), 0, 1);
Y(:,2) = params(1)*normcdf(params(2)*X(:,2)+params(4), 0, 1);
% 
% Y(:,1) = params(1)./(1+exp(params(2)*(X(:,1)+params(3))));
% Y(:,2) = params(1)./(1+exp(params(2)*(X(:,2)+params(4))));
end