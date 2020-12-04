function params = fitNonlinearityParams(xBin, yBin)
% params = fitNonlinearityParams(xBin, yBin)

opts = statset('nlinfit');
opts.MaxFunEvals = 1e5;
opts.MaxIter = 1e5;

% W = 0.2*ones(size(yBin)); W(end-4:end) = 1;

% params = nlinfit(xBin,yBin,@models.ln.outputNonlinearity,[max(yBin)*3 0.2 -1.5 min(yBin)],opts,'Weights',W);
params = nlinfit(xBin,yBin,@outputNonlinearity,[max(yBin)*3 0.2 -1.5 min(yBin)],opts);
