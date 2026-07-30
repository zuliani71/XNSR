function weights = KonnoOhmachiWeights(f,fc,b)
%KONNOOHMACHIWEIGHTS Build normalized reusable Konno-Ohmachi weights.
%   W = KonnoOhmachiWeights(F,FC,B) returns a numel(F)-by-numel(FC)
%   matrix. Every column contains the normalized smoothing window for one
%   center frequency. A filtered spectrum is evaluated as W.' * X.

f = f(:);
fc = fc(:);

weights = KonnoOhmachiSmoothingWindow(f,fc,b);
normalization = sum(weights,1);

if any(~isfinite(normalization) | normalization <= 0)
    error('XNSR:KonnoOhmachiWeights:InvalidNormalization', ...
        'Every center frequency must have a finite, positive window sum.');
end

weights = weights ./ normalization;
end
