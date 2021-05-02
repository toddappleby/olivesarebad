function Y = multiVarNL(params, X)

%2 and 4 are gain
%3 and 5 are hz

Y(:,1) = params(1)*normcdf(params(2)*X(:,1)+params(3), 0, 1);
Y(:,2) = params(1)*normcdf(params(4)*X(:,2)+params(5), 0, 1);
% 
% Y(:,1) = params(1)./(1+exp(params(2)*(X(:,1)+params(3))));
% Y(:,2) = params(1)./(1+exp(params(4)*(X(:,2)+params(5))));
end