function [errCatch, paramsCatch] = plotNonlinearityComparison(typeName, mrComparison, varargin)
% function plotNonlinearityComparison(xBinLo, yBinLo, xBinHi, yBinHi, varargin)
% plotNonlinearityComparison(xBinLo, yBinLo, xBinHi, yBinHi)
% Fits and plots nonlinearity comparisons with overlays and scaling
%
% Parameters:
% xBinLo, yBinLo: x/y values for low condition (red curve)
% xBinHi, yBinHi: x/y values for high condition (blue curve)
% Optional parameters:
% 'adaptationType': 'gain','horizontal','gain+horizontal' (default)
% 'method': 'function','extrapolation' (default),'multifit'

    % Parse inputs
    p = inputParser();
    p.addParameter('adaptationType', 'gain+horizontal', @(x)ischar(x));
    p.addParameter('method', 'function', @(x)ischar(x));
    p.parse(varargin{:});

    % open cell type folder
    exportFolder = getDirectory();
    protocolToAnalyze = 'Motion And Noise_Mike';
    typesFolder = strcat(exportFolder,protocolToAnalyze,'/','celltypes/',typeName);
    cd(typesFolder)
    typeDir = dir;
    fileSet = char(string({typeDir.name}));
    fileIndex = fileSet(:,:,fileSet(:,1,:) == typeName(1));
 
    %loop through each cell, plot results, grab errors
    
    figure('Position', [-1350 -850 800 1600]); %for vertical 2nd monitor
    
    grabExample = true;
    exampleNumber = 4;
    
% figure('Position', [100 100 800 600])
    for d = 1:size(fileIndex,3)
    
   

        load(strtrim(fileIndex(:,:,d)))

        % check nonlinearity sequence (most are 1. random 2. motion 3. static)
        if strcmp(analysisParams.classLabels(1),'random')  
            xRandom = analysisParams.xBin(1,:);
            yRandom = analysisParams.yBin(1,:);
            xMotion = analysisParams.xBin(2,:);
            yMotion = analysisParams.yBin(2,:);
            xStatic = analysisParams.xBin(3,:);
            yStatic = analysisParams.yBin(3,:);
        else
            xRandom = analysisParams.xBin(2,:);
            yRandom = analysisParams.yBin(2,:);
            xMotion = analysisParams.xBin(1,:);
            yMotion = analysisParams.yBin(1,:);
            xStatic = analysisParams.xBin(3,:);
            yStatic = analysisParams.yBin(3,:);
        end

        % Fit the nonlinearities
        if mrComparison 
        [nlParams, err, r2] = fitNonlinearityGainChange(xMotion, yMotion, xRandom, yRandom, ...
            'adaptationType', p.Results.adaptationType, ...
            'method', p.Results.method);
        else
        [nlParamsMS, errMS, r2MS] = fitNonlinearityGainChange(xMotion, yMotion, xStatic, yStatic, ...
            'adaptationType', p.Results.adaptationType, ...
            'method', p.Results.method);
        
         [nlParamsRS, errRS, r2RS] = fitNonlinearityGainChange(xRandom, yRandom, xStatic, yStatic, ...
            'adaptationType', p.Results.adaptationType, ...
            'method', p.Results.method);
        end

        % Create figure with subplots
        
         
        % Plot 1: Original curves overlay
        subplot(size(fileIndex,3),2,1+(2*(d-1)));
       
%         subplot(1,2,1+(2*(1-1)));
        plot(xRandom, yRandom, 'b-', 'LineWidth', 2);
        hold on;
        plot(xMotion, yMotion, 'r-', 'LineWidth', 2);
        title('Overlay of Original Nonlinearities');
        xlabel('Filtered Stim Input');
        ylabel('Response (Hz)');
        grid on;

        % Plot 2: Scaled comparison
        subplot(size(fileIndex,3),2,2+(2*(d-1)));
        
%         subplot(1,2,2+(2*(1-1)));
        plot(xRandom, yRandom, 'b-', 'LineWidth', 2);
        hold on;

        % Apply the fitted transformation based on adaptation type
        switch p.Results.adaptationType
            case 'gain'
                xScaled = xMotion/nlParams(1);
                yScaled = yMotion;
            case 'horizontal'
                xScaled = xMotion - nlParams;
                yScaled = yMotion;
            case 'gain+horizontal'
                xScaled = xMotion/nlParams(1) + nlParams(2);
                yScaled = yMotion;
            case 'gain+vertical'
                xScaled = xMotion/nlParams(1);
                yScaled = yMotion - nlParams(2);
        end

        % Plot scaled curve
        plot(xScaled, yScaled, 'r--', 'LineWidth', 2);
        title('Scaled Comparison',p.Results.adaptationType);
        xlabel('Filtered Input');
        ylabel('Response');
        grid on;
        sgtitle(strcat('Adaptation Type: ',p.Results.adaptationType));

%         Add legend
        legend('Random Condition', 'Motion (Scaled) Condition', 'Location', 'northwest');

%        %Display fit quality metrics
%         if ~isempty(err)
%             fprintf('Fit Error (MSE): %.4f\n', err);
%         end
%         if ~isempty(r2)
%             fprintf('R-squared: %.4f\n', r2);
%         end
%         
        errCatch(d) = err;
        paramsCatch(d) = nlParams;
    end
    
end