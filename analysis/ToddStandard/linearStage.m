function x = linearStage(tparams, stimulus, sampleRate)
% x = linearStage(tparams, stimulus, sampleRate)

% Get the dimensionality of the stimulus.
d = ndims(stimulus);

% Make sure the stimulus is a row.
if d == 2 && size(stimulus,2) == 1
    stimulus = stimulus';
end

% Pad the front and back with zeros to avoid wrap effects.
padLength = floor(0.5 * sampleRate);

% Create the temporal filter.
t = (0 : size(stimulus,d)-1)/sampleRate;
tfilter = models.ln.linearFilterFunction(tparams, t);
% Scale to unit vector.
tfilter = tfilter / norm(tfilter);
tfilter = tfilter / sum(abs(tfilter));


if d == 2
    x = zeros(size(stimulus) + [0 padLength*2]);
    for k = 1 : size(stimulus,1)
        % Pad and convolve.
        x(k,:) = real(ifft( fft([tfilter,zeros(1,padLength*2)]) .* fft([zeros(1,padLength),stimulus(k,:),zeros(1,padLength)]) ));
    end
    % Remove the padding.
    x = x(:,padLength+(1 : length(stimulus)));
else
    x = zeros(size(stimulus) + [0 0 padLength*2]);
    for k = 1 : size(stimulus,1)
        for m = 1 : size(stimulus,2)
            % Pad and convolve.
            x(k,m,:) = real(ifft( fft([tfilter,zeros(1,padLength*2)]) .* fft([zeros(1,padLength),squeeze(stimulus(k,m,:))',zeros(1,padLength)]) ));
        end
    end
    % Remove the padding.
    x = x(:,:,padLength+(1 : length(stimulus)));
end