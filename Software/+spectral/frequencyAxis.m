function f = frequencyAxis(n, fs)
%FREQUENCYAXIS Non-negative FFT frequencies for N samples at rate FS.

arguments
    n (1,1) double {mustBeInteger, mustBePositive}
    fs (1,1) double {mustBeFinite, mustBePositive}
end

f = (0:floor(n / 2)).' * (fs / n);
end
