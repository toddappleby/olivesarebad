function [nlParams,err,r2] = fitNonlinearityGainChange(xBinLo,yBinLo,xBinHi,yBinHi,varargin)
% [nlParams,err,r2] = fitNonlinearityGainChange(xBinLo,yBinLo,xBinHi,yBinHi,varargin)
%
% Parameters:
% xBinLo/yBinLo = x/y values for your nonlinearity for one condition
% xBinHi/yBinHi = x/y values for your nonlinearity for second condition
% 'adaptationType' = 'gain', 'horizontal', 'gain+horizontal', 'gain+vertical', 'gain+horizontal+vertical'
% 'method':
%       'function' = parameteric (cumulative Normal) fit based on Chichilnisky et al., 2001
%       'extrapolation' = non-parameteric (extrapolation) fit
%       'multifit' = performs both 'function' and 'extrapolation' fits and selects the one with the lowest error.
% 
% Returns:
% nlParams returns the x-axis scale or offset, etc, 
% err is the MSE, and r2 is the r^2 value for the fit.
%
% Examples:
% [nlParams,err,r2] = fitNonlinearityGainChange(xBinLo, yBinLo, xBinHi, yBinHi, 'adaptationType', 'gain', 'method', 'extrapolation')
% [nlParams,err,r2] = fitNonlinearityGainChange(xBinLo, yBinLo, xBinHi, yBinHi, 'adaptationType', 'horizontal', 'method', 'extrapolation')

    ip = inputParser();
    ip.addParameter('adaptationType', 'gain+horizontal', @(x)ischar(x));
    ip.addParameter('method', 'function', @(x)ischar(x));
    ip.parse(varargin{:});

    adaptationType = ip.Results.adaptationType;
    method = ip.Results.method;

    opts = statset('nlinfit');
    opts.MaxFunEvals = 1e4; 
    opts.MaxIter = 1e4;
    r2 = [];

    switch method
        case 'function'
            [nlParams,err,r2] = fitGainFunctionMethod(xBinLo,yBinLo,xBinHi,yBinHi,adaptationType,opts);
        case 'multifit'
            [nlParams,err] = multifitExtrapMethod(xBinLo,yBinLo,xBinHi,yBinHi,adaptationType);
        otherwise
            [nlParams,err] = fitGainExtrapMethod(xBinLo,yBinLo,xBinHi,yBinHi,adaptationType,opts);
    end
end

function [nlParams, err, r2] = fitGainFunctionMethod(xBinLo,yBinLo,xBinHi,yBinHi,adaptationType,opts)
    nlfun = @(a,x)(a(1)*normcdf(a(2)*x+a(3),0,1)+a(4));

    hiParams = nlinfit(xBinHi(:),yBinHi(:),nlfun,[70 2 -1.5 -5],opts);

    switch adaptationType
        case 'gain+vertical' % Fit gain change and vertical offset between high and low contrast.
            nlfun2 = @(a,x)(hiParams(1)*normcdf(a(1)*hiParams(2)*x+hiParams(3),0,1)+hiParams(4)+a(2));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,nlfun2,[1 0],opts);
        case 'gain+horizontal' % Gain change and horizontal shift
            nlfun2 = @(a,x)(hiParams(1)*normcdf(a(1)*hiParams(2)*x+(a(2)+hiParams(3)),0,1)+hiParams(4));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,nlfun2,[1 0],opts);
        case 'gain+horizontal+vertical' % Gain change plus horizontal and vertical shifts
            nlfun2 = @(a,x)(interp1(xBinHi/a(1)+a(2),yBinHi,x,'linear','extrap') + a(3));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,nlfun2,[0.7 0 0],opts);
        case 'horizontal' % Horizontal shift only.
            nlfun2 = @(a,x)(hiParams(1)*normcdf(hiParams(2)*x+(a+hiParams(3)),0,1)+hiParams(4));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,nlfun2,1,opts);
        otherwise % Gain change only
%             nlfun2 = @(a,x)(hiParams(1)*normcdf(a(1)*hiParams(2)*x+(hiParams(3)),0,1)+hiParams(4));
            nlfun2 = @(a,x)(interp1(xBinHi/a,yBinHi,x,'linear','extrap'));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,nlfun2,0.7,opts);
    end

    % Get the r^2 value.
    
    r = corrcoef([yBinHi(:)',yBinLo(:)'], [nlfun(hiParams,xBinHi),nlfun2(nlParams,xBinLo)]);
    
    r2 = r(2);
end

function [nlParams,err] = fitGainExtrapMethod(xBinLo,yBinLo,xBinHi,yBinHi,adaptationType,opts)
    % Get the initial guess for gain change.
%     g = yBinHi(:) \ yBinLo(:);
    
    g = 0.6;
    
    
    switch adaptationType
        case 'gain+vertical' % Fit gain change and vertical offset between high and low contrast.
            modelfun = @(a,x)(interp1(xBinHi*a(1),yBinHi,x,'linear','extrap') + a(2));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,modelfun,[g 0],opts);
        case 'gain+horizontal' % Gain change and horizontal shift
            
            x = xBinHi(yBinHi <= max(yBinLo));
            y = yBinHi(yBinHi <= max(yBinLo));
            % Project the y-axis onto the x-axis to simplify the
            % computation.
            xNew = linspace(min([x(:); xBinLo(:)]),max([x(:); xBinLo(:)]),100);
            yIntH = interp1(x, y,xNew,'linear','extrap');
            yIntL = interp1(xBinLo, yBinLo,xNew,'linear','extrap');
            xOffset = xNew(find(abs(yIntH-yIntL) == min(abs(yIntH-yIntL)),1));
            
            pL = polyfit(yIntL(xNew >= xOffset),xNew(xNew >= xOffset),1);
            pH = polyfit(yIntH(xNew >= xOffset),xNew(xNew >= xOffset),1);
            
            
            modelfun = @(a,x)(interp1(pL(1)/pH(1)*xBinHi+a(1),yBinHi,x,'linear','extrap'));
            [tmpParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,modelfun,0,opts);
            nlParams = [pL(1)/pH(1), tmpParams(1)];
        case 'gain+horizontal+vertical' % Gain change plus horizontal and vertical shifts
            modelfun = @(a,x)(interp1(a(1)*xBinHi+a(2),yBinHi,x,'linear','extrap') + a(3));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,modelfun,[g 0.02 0],opts);
        case 'horizontal' % Horizontal shift only.
            modelfun = @(a,x)(interp1(xBinHi+a(1),yBinHi,x,'linear','extrap'));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,modelfun,0,opts);
        otherwise % Gain change only
            modelfun = @(a,x)(interp1(a*xBinHi,yBinHi,x,'linear','extrap'));
            [nlParams,~,~,~,err] = nlinfit(xBinLo,yBinLo,modelfun,g,opts);
    end
    
    
   
    
end

% Method for looping through each of the possible fitting methods.
function [nlParams,err] = multifitExtrapMethod(xBinLo,yBinLo,xBinHi,yBinHi,adaptationType)
    % Get the initial guess for gain change.
    g = yBinHi(:) \ yBinLo(:);
    
    options = optimoptions('lsqcurvefit',...
        'Display','off',...
        'MaxFunctionEvaluations',1e5,...
        'MaxIterations',1e5);
    
    nlP = zeros(2,2);
    err = zeros(2,1);
    for k = 1 : 2
        switch k
            case 1 % Horizontal shift only.
                if strcmpi(adaptationType,'horizontal')
                    lb=-0.5; ub=0.5;
                    modelfun = @(a,x)(interp1(xBinHi+a(1),yBinHi,x,'linear','extrap'));
                    [nl,err(k)] = lsqcurvefit(modelfun,0,xBinLo,yBinLo,lb,ub,options);
                    nlP(k,2) = nl;
    %                 [nl,~,~,~,err(k)] = nlinfit(xBinLo,yBinLo,modelfun,0,opts);
                else % Gain change only
                    lb=0.5; ub=1.5;
                    modelfun = @(a,x)(interp1(a*xBinHi,yBinHi,x,'linear','extrap'));
                    [nlP(k,1),err(k)] = lsqcurvefit(modelfun,g,xBinLo,yBinLo,lb,ub,options);
    %                 [nlParams(k,1),~,~,~,err(k)] = nlinfit(xBinLo,yBinLo,modelfun,g,opts);
                end
            case 2 % Gain change and horizontal shift
                lb=[0.5 -0.5]; ub=[1.5 0.5];
                modelfun = @(a,x)(interp1(a(1)*xBinHi+a(2),yBinHi,x,'linear','extrap'));
                if strcmpi(adaptationType,'horizontal')
                    [nlP(k,:),err(k)] = lsqcurvefit(modelfun,[1 nlP(1,2)],xBinLo,yBinLo,lb,ub,options);
                else
                    [nlP(k,:),err(k)] = lsqcurvefit(modelfun,[nlP(1,1) 0],xBinLo,yBinLo,lb,ub,options);
                end
%                 [nlParams(k,:),~,~,~,err(k)] = nlinfit(xBinLo,yBinLo,modelfun,[nlParams(2,1) nlParams(1,2)],opts);
        end
    end
    err(2) = err(2)*2;
    err = err / length(xBinLo);
    
    index = find(err == min(err),1);
    nlParams = nlP(index,:);
   
    err = err(index);
end
