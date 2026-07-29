function H = halfSpectrum(X, dim)
%HALFSPECTRUM Extract the non-negative-frequency half of an FFT.
%   H = spectral.halfSpectrum(X) operates along the first dimension.
%   H = spectral.halfSpectrum(X, DIM) operates along dimension DIM.
%
%   No scaling or amplitude normalization is applied.

arguments
    X {mustBeNumeric}
    dim (1,1) double {mustBeInteger, mustBePositive} = 1
end

n = size(X, dim);
indices = repmat({':'}, 1, max(ndims(X), dim));
indices{dim} = 1:(floor(n / 2) + 1);
H = X(indices{:});
end
