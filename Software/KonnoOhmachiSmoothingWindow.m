function WIN = KonnoOhmachiSmoothingWindow(f,fc,b)
%KONNOOHMACHISMOOTHINGWINDOW Build Konno-Ohmachi smoothing windows.
%   WIN = KONNOOHMACHISMOOTHINGWINDOW(F,FC,B) evaluates the windows on
%   frequency vector F for each center frequency in FC. A scalar FC
%   produces a column vector; a vector FC produces one column per center.
%   B controls the logarithmic bandwidth.
%
%   Example:
%       f = 0.1:0.1:100;
%       win = KonnoOhmachiSmoothingWindow(f,1:5,40);
%       semilogx(f,win);
%
%   Originally written by D. Zuliani.

% Normalize frequency inputs to column vectors.
f=f(:);
fc=fc(:);
% Expand frequency and center-frequency grids.
[fc,f]  =   meshgrid(fc,f);
% Evaluate the Konno-Ohmachi argument.
X       = b*log10(f./fc);
SINX    = sin(X);
% Define sin(x)/x as one at x = 0.
I = find(SINX==0 & X==0);
X(I)    =   1;
SINX(I) =   1;
% Frequencies at zero produce -Inf in log10 and receive zero weight.
I = find(isinf(X));
SINX(I) =   0;
% Apply the fourth-power Konno-Ohmachi window.
WIN     = (SINX./X).^4;
