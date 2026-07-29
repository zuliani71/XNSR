function xout = fft2ft(xin, varargin)
%FFT2FT Legacy-compatible FFT to one-sided-spectrum conversion.
%   H = FFT2FT(X) extracts the raw non-negative-frequency half.
%   OUT = FFT2FT(X,FS) returns a one-sided amplitude spectrum using the
%   historical vector or page-based output layout.
%
%   This XNSR compatibility entry point delegates the numerical work to
%   the vendored +spectral package distributed in Software/+spectral.

narginchk(1, 2);
if isvector(xin)
    % Historical FFT2FT behavior: vectors are always processed and
    % returned as columns, irrespective of their input orientation.
    xin = xin(:);
end

if nargin == 1
    xout = spectral.halfSpectrum(xin, 1);
    return
end

fs = varargin{1};
[amplitude, frequency] = spectral.oneSidedAmplitudeFromFFT(xin, fs, 1);

if isvector(xin)
    xout = [frequency, amplitude(:)];
    return
end

frequencyShape = ones(1, max(ndims(amplitude), 2));
frequencyShape(1) = numel(frequency);
frequency = reshape(frequency, frequencyShape);
frequency = repmat(frequency, [1, size(amplitude, 2), ...
    ones(1, max(ndims(amplitude) - 2, 0))]);
xout = cat(3, frequency, amplitude);
end
