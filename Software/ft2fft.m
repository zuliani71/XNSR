function xout = ft2fft(xin, varargin)
%FT2FFT Legacy-compatible reconstruction of a full FFT spectrum.
%   X = FT2FFT(H) assumes Hermitian symmetry and odd original length.
%   X = FT2FFT(H,SYMMETRY) also assumes odd original length.
%   X = FT2FFT(H,SYMMETRY,PARITY) selects even ('e') or odd ('o') N.
%
%   SYMMETRY is 'o' for Hermitian (conjugate) or 'e' for plain mirror.
%   Matrices and N-D arrays are reconstructed along their first dimension.

narginchk(1, 3);
symmetryFlag = 'o';
parityFlag = 'o';
if nargin >= 2
    symmetryFlag = varargin{1};
end
if nargin == 3
    parityFlag = varargin{2};
end

symmetryFlag = validatestring(symmetryFlag, {'o','e'}, mfilename, ...
    'symmetry flag');
parityFlag = validatestring(parityFlag, {'o','e'}, mfilename, ...
    'parity flag');

halfLength = size(xin, 1);
if parityFlag == 'e'
    n = 2 * max(halfLength - 1, 0);
else
    n = 2 * halfLength - 1;
end
if n < 1
    error("spectral:ft2fft:EmptyInput", "The half spectrum cannot be empty.");
end

if symmetryFlag == 'o'
    symmetry = "hermitian";
else
    symmetry = "mirror";
end
xout = spectral.rebuildSpectrum(xin, n, symmetry, 1);
end
