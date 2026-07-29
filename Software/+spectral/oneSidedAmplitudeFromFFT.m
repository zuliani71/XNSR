function [A, f] = oneSidedAmplitudeFromFFT(X, fs, dim)
%ONESIDEDAMPLITUDEFROMFFT One-sided amplitude spectrum from FFT values.
%   [A,F] = spectral.oneSidedAmplitudeFromFFT(X,FS) operates along the
%   first dimension. A sinusoid on an FFT bin retains its peak amplitude.
%   DC and, for even-length inputs, Nyquist are not doubled.

arguments
    X {mustBeNumeric}
    fs (1,1) double {mustBeFinite, mustBePositive}
    dim (1,1) double {mustBeInteger, mustBePositive} = 1
end

n = size(X, dim);
A = spectral.halfSpectrum(X, dim) / n;

indices = repmat({':'}, 1, max(ndims(A), dim));
if mod(n, 2) == 0
    interior = 2:(size(A, dim) - 1);
else
    interior = 2:size(A, dim);
end
if ~isempty(interior)
    indices{dim} = interior;
    A(indices{:}) = 2 * A(indices{:});
end

f = spectral.frequencyAxis(n, fs);
end
