function analyze_gain_changes()
    % Define stimulus conditions
    mean_low = 100;  % Example low mean intensity
    mean_high = 200; % Example high mean intensity
    contrast_low = 0.1;  % Example low contrast
    contrast_high = 0.2; % Example high contrast

    % Generate x values (representing filtered input)
    x = linspace(-5, 5, 1000);

    % Define two nonlinear curves representing RGC responses
    y_low = 1 ./ (1 + exp(-(x + 2)));  % Response at low mean/contrast
    y_high = 1 ./ (1 + exp(-0.7 * (x - 1)));  % Response at high mean/contrast (reduced gain)

    % Optimize scaling (gain change only)
    options = optimset('Display', 'iter');
    [alpha, min_error] = fminsearch(@(a) objective_function(x, y_low, x, y_high, a), 1, options);

    % Apply the optimal scaling to y_low
    y_low_scaled = interp1(x*alpha, y_low, x, 'linear', 'extrap');

    % Calculate gamma indices
    gamma_mean = calculate_gamma(alpha, mean_high, mean_low);
    gamma_contrast = calculate_gamma(alpha, contrast_high, contrast_low);

    % Create plots
    plot_results(x, y_low, y_high, x*alpha, y_low_scaled, alpha);

    % Display results
    display_results(alpha, min_error, gamma_mean, gamma_contrast, mean_high, mean_low, contrast_high, contrast_low);
end

function error = objective_function(x, y1, x2, y2, alpha)
    y1_scaled = interp1(x*alpha, y1, x, 'linear', 'extrap');
    
    % Calculate error only in the overlapping range
    valid = ~isnan(y1_scaled) & ~isnan(y2);
    error = sum((y1_scaled(valid) - y2(valid)).^2);
end

function gamma = calculate_gamma(alpha, high, low)
    gamma = (alpha - 1) / (high/low - 1);
end

function plot_results(x, y_low, y_high, x_scaled, y_low_scaled, alpha)
    figure;

    % Original curves
    subplot(2,2,1);
    plot(x, y_low, 'b-', x, y_high, 'r-', 'LineWidth', 2);
    title('Original Nonlinearities');
    xlabel('Filtered Input'); ylabel('Spike Rate');
    legend('Low Condition', 'High Condition');

    % Overlay
    subplot(2,2,2);
    plot(x, y_low, 'b-', x, y_high, 'r-', 'LineWidth', 2);
    title('Overlay');
    xlabel('Filtered Input'); ylabel('Spike Rate');
    legend('Low Condition', 'High Condition');

    % Scaled to match
    subplot(2,2,3);
    plot(x, y_high, 'r-', x_scaled, y_low, 'b-', 'LineWidth', 2);
    title('Scaled to Match');
    xlabel('\alpha * Filtered Input'); ylabel('Spike Rate');
    legend('High Condition', 'Low Condition (scaled)');

    % Scaled overlay
    subplot(2,2,4);
    plot(x, y_high, 'r-', x, y_low_scaled, 'b-', 'LineWidth', 2);
    title('Scaled Overlay');
    xlabel('Filtered Input'); ylabel('Spike Rate');
    legend('High Condition', 'Low Condition (scaled)');
end

function display_results(alpha, min_error, gamma_mean, gamma_contrast, mean_high, mean_low, contrast_high, contrast_low)
    fprintf('Optimal scaling factor (α): %.4f\n', alpha);
    fprintf('Minimum squared error: %.6f\n', min_error);
    fprintf('γ_mean: %.4f\n', gamma_mean);
    fprintf('γ_contrast: %.4f\n', gamma_contrast);
    fprintf('Fractional change in mean: %.4f\n', mean_high/mean_low - 1);
    fprintf('Fractional change in contrast: %.4f\n', contrast_high/contrast_low - 1);
    fprintf('Fractional change in gain: %.4f\n', alpha - 1);
end