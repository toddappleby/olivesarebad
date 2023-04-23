function Y = multiGainNL(params, X)




Y(:,1) = params(1)*normcdf(params(3)*X(:,1)+params(2), 0, 1);
Y(:,2) = params(1)*normcdf(params(4)*X(:,2)+params(2), 0, 1);
Y(:,3) = params(1)*normcdf(params(5)*X(:,2)+params(2), 0, 1);


% Y(:,1) = params(1)./(1+exp(params(3)*(X(:,1)+params(2))));
% Y(:,2) = params(1)./(1+exp(params(4)*(X(:,2)+params(2))));
end