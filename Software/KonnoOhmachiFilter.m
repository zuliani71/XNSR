function [Y] = KonnoOhmachiFilter(X,f,fc,b,varargin)
%
% Y = KonnoOhmachiFilter(X,f,fc,b,varargin)
% filters the input spectra X (both matrix or vector)
% using the KonnoOhmachi method.
% X can be:
% - a column vector that includes a signal single
%   left side FFT;
% - a 2D matrix whose columns are signal single left
%   side FFTs;
% - a 3D matrix whose pages accommodate 2D matrix as
%    formatted as before.
% f is the frequency vector (usually a column) where the
% FFTS have been evaluated.
% fc is the center frequency for the KonnoOhmachi filter
% evaluation. fc can be both a scalar or a vector.
% b is the Konno Omachi coefficient strcitly related to
% the bandwidth of the window.
% Konno and Ohmachi smoothing function is strongly
% recommended for frequency analysis because it ensures
% a constant number points at low and high frequency.
%
% e.g. 
%      f = (0.1:0.1:10)'; X = f.^2; fc = (1:5)';  
%      Y = KonnoOhmachiFilter(X,f,fc,40);
%      semilogx(fc,Y);
% e.g. 
%      f = (0.1:0.1:10)'; X = [f,f.^2]; fc = (1:5)';  
%      [Y] = KonnoOhmachiFilter(X,f,fc,40);
%      semilogx(f,X); hold on; semilogx(fc,Y,'r');
%
% Made by D. Zuliani 2013/08/20

%
% MATRIX DIMENTION CHECK
NDIMS = ndims(X);
switch NDIMS
    case 2
        % NOTHING TO DO
    case 3
        % ACCOMODATING THE 3D STRUCTURE IN
        % A 2D MATRIX
        XSIZE = size(X);
        X = reshape(X,[XSIZE(1) XSIZE(2)*XSIZE(3) 1]);
    otherwise
end
%
% WORKING WITH COLUMN VECTORS
fc  = fc(:);
f   = f(:);
%
% GETTING CRUCIAL DIMs
M   = length(fc);   % number of KonnoHomachi center frequencies
N   = length(f);    % number of input frequencies
K   = size(X,2);    % number of input signals to filter
%
% CREATING THE KONNOOHMACHI COEFFICIENTS.
% They are accommodated inside a NxM matrix: N=lenght(f), M=length(fc).
WIN = KonnoOhmachiSmoothingWindow(f,fc,b);
%
% CREATING A SINGLE COLUMN OF KONNO COEFFS.
% All of them are clustered in groups fc1, fc2, ..., fcM. Each
% group is N length and includes  the KonnoOhmachi coeffs. for
% a given central frequency fc
WIN = reshape(WIN,size(WIN,2)*size(WIN,1),1);
%
% REPEATING INPUT SIGNAL MATRIX ALONG M ROWS
% Each matrix is made by the splitted signals Sig1,Sig2,...,SigK
% vectors.
SIG = repmat(X,M,1);
%
% ACCOMODATING THE COLUMNS OF X IN A SINGLE COLUMN
% The signal column vectors of the SIG matrix are
% reshaped in a single column made by M*K*N rows.
% the 1st signal vector is made by the 1st N values,
% it is repeated M times and this is the 1st cluster
% including multiple copies of the 1st vector. The
% 2nd cluster is made in the same way for the 2nd
% signal.
SIG = reshape(SIG,M*K*N,1);
%
% REPEATING THE KONNOMACHI COEFF. VECTOR.
% the WIN vector is repeated to match the exact number
% of SIG rows. The fc1,fc2,...,fcM groups are repeated
% K times along rows.
WIN = repmat(WIN,K,1);
%
% BUILDING THE SIMPLE SIG*WIN PRODUCT
% This is the product element by element between the
% signal vector and the column made by the fc groups.
% The operation is the 1st part of the KonnoOhmachi
% filtering;
Y = SIG.*WIN;
%
% WORKING WITH SUBINDEX
% They are needed to implement the full KonnoOhmachi
% method using accumarray
SUBS = 1:1:M*K;
SUBS = repmat(SUBS,N,1);
SUBS = reshape(SUBS,N*M*K,1);
%
% NORMALIZATION FACTOR
NFACT = accumarray(SUBS,WIN);
%
% KONNOMACHI FILTERING LAST PART.
% The element by element products between the signal
% vector and the column made by the fc groups are finally
% summed group by group following the rule imposed by SUBS.
% The normalization factor is applied too.
Y = accumarray(SUBS,Y)./NFACT;
%
% RESHAPING
% Y is a matrix made by K column vectors of M size.
Y = reshape(Y,M,K);
%
% REBUILD A 3D MATRIX IF NEEDED
switch NDIMS
    case 2
        % NOTHING TO DO
    case 3
        % ACCOMODATING THE 2D STRUCTURE IN
        % A 3D MATRIX
        Y = reshape(Y,[M XSIZE(2) XSIZE(3)]);
    otherwise
end