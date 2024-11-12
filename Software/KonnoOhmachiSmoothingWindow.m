function WIN = KonnoOhmachiSmoothingWindow(f,fc,b)
% 
% [WIN] = KonnoOhmachiSmoothingWindow(f,fc,b,varargin)
% calculates the Konno Omachi smoothing window WIN
% around the given center frequency fc respect to the
% frequency vector f.
% fc can be both a scalar or a vector.
% If fc is a scalar, WIN is a column vector the same size of f,
% if fc is a vector, WIN is a matrix with a number of rows
% equal to the length of f and a number of columns equal to
% the length of fc.
% b is the Konno Omachi coefficient strcitly related to
% the band width of the window.
% Konno and Ohmachi smoothing function is strongly
% recommended for frequency analysis because it ensures
% a constant number points at low and high frequency.
%
% e.g. WIN = KonnoOhmachiSmoothingWindow(0.1:0.1:100,1:5,40);
%      semilogx(0.1:0.1:100,WIN);
%
% Made by D. Zuliani 2013/08/20

% 
% Build a column frequency vector
f=f(:);
%
% Build a column central frequency vector
fc=fc(:);
%
% build the f and fc matrix
[fc,f]  =   meshgrid(fc,f);
%
% Construct the main pieces for the the Konno Omachi rule
X       = b*log10(f./fc);
SINX    = sin(X);
%
% function sin(X)/X=1 for X=0
I = find(SINX==0 & X==0);
X(I)    =   1;
SINX(I) =   1;
%
% WARNING if f = 0 then log10 is -Inf
I = find(isinf(X));
SINX(I) =   0;
%
% Yeld the Konno Omachi window
WIN     = (SINX./X).^4;