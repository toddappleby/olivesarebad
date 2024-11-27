function plotNonlinearityComparisonSimple(xBinMotion, yBinMotion, xBinRandom, yBinRandom, varargin)
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

    
    % Fit the nonlinearities
    [nlParams, err, r2] = fitNonlinearityGainChange(xBinMotion, yBinMotion, xBinRandom, yBinRandom, ...
        'adaptationType', p.Results.adaptationType, ...
        'method', p.Results.method);
    
    % Create figure with subplots
    figure('Position', [100 100 800 400]);
    
    % Plot 1: Original curves overlay
    subplot(1,2,1);
    plot(xBinRandom, yBinRandom, 'b-', 'LineWidth', 2);
    hold on;
    plot(xBinMotion, yBinMotion, 'r-', 'LineWidth', 2);
    title('Overlay of Original Nonlinearities');
    xlabel('Filtered Stim Input');
    ylabel('Response (Hz)');
    grid on;
    
    % Plot 2: Scaled comparison
    subplot(1,2,2);
    plot(xBinRandom, yBinRandom, 'b-', 'LineWidth', 2);
    hold on;
    
    % Apply the fitted transformation based on adaptation type
    switch p.Results.adaptationType
        case 'gain'
            xScaled = xBinMotion/nlParams(1);
            yScaled = yBinMotion;
        case 'horizontal'
            xScaled = xBinMotion - nlParams;
            yScaled = yBinMotion;
        case 'gain+horizontal'
            xScaled = xBinMotion/nlParams(1) + nlParams(2);
            yScaled = yBinMotion;
        case 'gain+vertical'
            xScaled = xBinMotion/nlParams(1);
            yScaled = yBinMotion - nlParams(2);
    end
    
    % Plot scaled curve
    plot(xScaled, yScaled, 'r--', 'LineWidth', 2);
    title('Scaled Comparison');
    xlabel('Filtered Input');
    ylabel('Response');
    grid on;
    
    % Add legend
    legend('Random Condition','Motion (Scaled) Condition', 'Location', 'northwest');
    
    % Display fit quality metrics
    if ~isempty(err)
        fprintf('Fit Error (MSE): %.4f\n', err);
    end
    if ~isempty(r2)
        fprintf('R-squared: %.4f\n', r2);
    end
end