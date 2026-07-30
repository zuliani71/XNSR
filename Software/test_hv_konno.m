% Regression tests for the standalone simplified H/V utility.

softwarePath = fileparts(mfilename('fullpath'));
addpath(softwarePath,'-begin');
rng(17);

sampleFrequency = 100;
konnoCoefficient = 40;
sampleCount = 1024;
baseSignal = 2 + randn(sampleCount,1);

% Identical components must yield H/V = 1 for even and odd FFT sizes.
for fftSize = [512,511]
    XYZ = [baseSignal,baseSignal,baseSignal];
    HV = hv_konno(XYZ,sampleFrequency,konnoCoefficient,fftSize);
    assert(size(HV,1) == floor(fftSize/2)+1);
    assert(size(HV,2) == 2);
    assert(all(isfinite(HV),'all'));
    assert(max(abs(HV(:,2)-1)) <= 100*eps);
    assert(abs(HV(end,1)-floor(fftSize/2)* ...
        sampleFrequency/fftSize) <= 10*eps(sampleFrequency));
end

% Verify the implementation against an explicit reference calculation.
XYZ = 1 + abs(randn(sampleCount,3));
fftSize = 512;
[amplitude,frequency] = spectral.oneSidedAmplitudeFromFFT( ...
    fft(XYZ,fftSize,1),sampleFrequency,1);
amplitude = abs(amplitude);
horizontal = sqrt(amplitude(:,1).*amplitude(:,2));
vertical = amplitude(:,3);
positive = frequency > 0;
weights = KonnoOhmachiWeights( ...
    frequency(positive),frequency(positive),konnoCoefficient);
reference = [horizontal,vertical];
reference(positive,:) = KonnoOhmachiFilter( ...
    reference(positive,:),frequency(positive),frequency(positive), ...
    konnoCoefficient,weights);
expected = [frequency,reference(:,1)./reference(:,2)];
actual = hv_konno(XYZ,sampleFrequency,konnoCoefficient,fftSize);
assert(max(abs(actual-expected),[],'all') <= 100*eps(max(expected,[],'all')));

% Negative time-domain samples must not create complex H/V values.
XYZ = randn(sampleCount,3);
XYZ(:,3) = XYZ(:,3) + 1;
actual = hv_konno(XYZ,sampleFrequency,konnoCoefficient,fftSize);
assert(isreal(actual));
assert(all(isfinite(actual),'all'));

fprintf('TUTTI I TEST HV_KONNO SONO SUPERATI\n');
