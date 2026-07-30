function y = triangFilter(x,f,fc,pcent,varargin)
%TRIANGFILTER Smooth spectra with variable-width triangular windows.
%   Y = TRIANGFILTER(X,F,FC,PERCENT) builds the smoothing weights and
%   applies them to vectors, matrices, or pages of spectra.
%
%   Y = TRIANGFILTER(X,F,FC,PERCENT,WEIGHTS) reuses a matrix returned by
%   TRIANGFILTERWEIGHTS. This form avoids rebuilding identical windows in
%   repeated calls, including calls made from a PARFOR loop.

narginchk(4,5);
validateattributes(x,{'numeric'},{'nonempty'},mfilename,'x',1);

inputSize = size(x);
inputDimensions = ndims(x);
assert(inputDimensions <= 3, ...
    'XNSR:TriangFilter:InvalidDimensions', ...
    'x must be a vector, matrix, or three-dimensional array.');
assert(size(x,1) == numel(f), ...
    'XNSR:TriangFilter:FrequencySizeMismatch', ...
    'The first dimension of x must match the number of input frequencies.');

if isempty(varargin)
    weights = triangFilterWeights(f,fc,pcent);
else
    weights = varargin{1};
    validateattributes(weights,{'numeric'}, ...
        {'2d','nonempty','finite'},mfilename,'weights',5);
    assert(isequal(size(weights),[numel(f),numel(fc)]), ...
        'XNSR:TriangFilter:WeightSizeMismatch', ...
        'weights must have size numel(f)-by-numel(fc).');
end

spectra = reshape(x,size(x,1),[]);
y = weights.' * spectra;

if inputDimensions == 3
    y = reshape(y,[numel(fc),inputSize(2),inputSize(3)]);
end
end
