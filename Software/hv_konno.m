function HV = hv_konno(XYZ,varargin)
%HV_KONNO Calculate a simplified Konno-Ohmachi H/V spectrum.
%   HV = HV_KONNO(XYZ) accepts an N-by-3 matrix containing the two
%   horizontal components followed by the vertical component. The output
%   columns contain frequency and smoothed H/V ratio.
%
%   HV = HV_KONNO(XYZ,FS,B,NFFT) selects sampling frequency FS
%   (default 100 Hz), Konno-Ohmachi coefficient B (default 40), and FFT
%   size NFFT (default 512).
%
%   This is a standalone simplified utility. XN_Cruncher performs its own
%   segmented, rotated and parallel H/V processing and does not call it.

narginchk(1,4);
validateattributes(XYZ,{'numeric'}, ...
    {'2d','ncols',3,'nonempty','finite'},mfilename,'XYZ',1);

sampleFrequency = 100;
konnoCoefficient = 40;
fftSize = 512;
if numel(varargin) >= 1 && ~isempty(varargin{1})
    sampleFrequency = varargin{1};
end
if numel(varargin) >= 2 && ~isempty(varargin{2})
    konnoCoefficient = varargin{2};
end
if numel(varargin) == 3 && ~isempty(varargin{3})
    fftSize = varargin{3};
end

validateattributes(sampleFrequency,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename,'FS',2);
validateattributes(konnoCoefficient,{'numeric'}, ...
    {'scalar','real','finite','positive'},mfilename,'B',3);
validateattributes(fftSize,{'numeric'}, ...
    {'scalar','integer','finite','positive'},mfilename,'NFFT',4);

% Build correctly normalized one-sided amplitude spectra for both
% horizontal components and the vertical component.
fullSpectrum = fft(XYZ,fftSize,1);
[amplitude,frequency] = spectral.oneSidedAmplitudeFromFFT( ...
    fullSpectrum,sampleFrequency,1);
amplitude = abs(amplitude);

horizontalAmplitude = sqrt(amplitude(:,1).*amplitude(:,2));
verticalAmplitude = amplitude(:,3);

% Konno-Ohmachi windows are defined for positive center frequencies.
% Preserve the unsmoothed DC value and smooth all strictly positive bins.
positive = frequency > 0;
smoothed = [horizontalAmplitude,verticalAmplitude];
if any(positive)
    weights = KonnoOhmachiWeights( ...
        frequency(positive),frequency(positive),konnoCoefficient);
    smoothed(positive,:) = KonnoOhmachiFilter( ...
        smoothed(positive,:),frequency(positive), ...
        frequency(positive),konnoCoefficient,weights);
end

HV = [frequency,smoothed(:,1)./smoothed(:,2)];
end
