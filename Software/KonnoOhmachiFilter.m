function Y = KonnoOhmachiFilter(X,f,fc,b,varargin)
%KONNOOHMACHIFILTER Smooth spectra using the Konno-Ohmachi method.
%   Y = KonnoOhmachiFilter(X,F,FC,B) filters spectra stored along the
%   first dimension of X. X can be a vector, matrix or N-D array.
%
%   Y = KonnoOhmachiFilter(X,F,FC,B,W) reuses the normalized weight matrix
%   W returned by KonnoOhmachiWeights. This form avoids rebuilding the same
%   window in repeated calls or inside PARFOR loops.
%
%   The output first dimension has numel(FC) elements. All remaining
%   dimensions are preserved.

narginchk(4,5);

f = f(:);
fc = fc(:);
n = numel(f);
m = numel(fc);

if size(X,1) ~= n
    error('XNSR:KonnoOhmachiFilter:FrequencySizeMismatch', ...
        'The first dimension of X (%d) must equal numel(f) (%d).', ...
        size(X,1),n);
end

if isempty(varargin)
    weights = KonnoOhmachiWeights(f,fc,b);
else
    weights = varargin{1};
    if ~isnumeric(weights) || ~isequal(size(weights),[n,m])
        error('XNSR:KonnoOhmachiFilter:InvalidWeights', ...
            'W must be a numeric %d-by-%d normalized weight matrix.',n,m);
    end
end

inputSize = size(X);
spectra = reshape(X,n,[]);

% Each output value is the normalized weighted sum of one input spectrum.
% Matrix multiplication avoids the large REPMAT/ACCUMARRAY temporaries used
% by the historical implementation.
filtered = weights.' * spectra;

outputSize = inputSize;
outputSize(1) = m;
Y = reshape(filtered,outputSize);
end
