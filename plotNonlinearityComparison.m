function plotNonlinearityComparison(xBinLo, yBinLo, xBinHi, yBinHi, varargin)
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
    [nlParams, err, r2] = fitNonlinearityGainChange(xBinLo, yBinLo, xBinHi, yBinHi, ...
        'adaptationType', p.Results.adaptationType, ...
        'method', p.Results.method);
    
    % Create figure with subplots
    figure('Position', [100 100 800 400]);
    
    % Plot 1: Original curves overlay
    subplot(1,2,1);
    plot(xBinHi, yBinHi, 'b-', 'LineWidth', 2);
    hold on;
    plot(xBinLo, yBinLo, 'r-', 'LineWidth', 2);
    title('Overlay of Original Nonlinearities');
    xlabel('Filtered Stim Input');
    ylabel('Response (Hz)');
    grid on;
    
    % Plot 2: Scaled comparison
    subplot(1,2,2);
    plot(xBinHi, yBinHi, 'b-', 'LineWidth', 2);
    hold on;
    
    % Apply the fitted transformation based on adaptation type
    switch p.Results.adaptationType
        case 'gain'
            xScaled = xBinLo/nlParams(1);
            yScaled = yBinLo;
        case 'horizontal'
            xScaled = xBinLo - nlParams;
            yScaled = yBinLo;
        case 'gain+horizontal'
            xScaled = xBinLo/nlParams(1) + nlParams(2);
            yScaled = yBinLo;
        case 'gain+vertical'
            xScaled = xBinLo/nlParams(1);
            yScaled = yBinLo - nlParams(2);
    end
    
    % Plot scaled curve
    plot(xScaled, yScaled, 'r--', 'LineWidth', 2);
    title('Scaled Comparison');
    xlabel('Filtered Input');
    ylabel('Response');
    grid on;
    
    % Add legend
    legend('High Condition', 'Low Condition', 'Location', 'northwest');
    
    % Display fit quality metrics
    if ~isempty(err)
        fprintf('Fit Error (MSE): %.4f\n', err);
    end
    if ~isempty(r2)
        fprintf('R-squared: %.4f\n', r2);
    end
end