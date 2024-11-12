function [y] = triangFilter(x,f,fc,pcent,varargin)
%
% y = triangFilter(x,f,fc,percent,varargin)
% filters the input spectra x (both matrix or vector)
% using a triangular mask of variable size. The size
% in frequency is a percent of each outupt frequency.
% x can be:
% - a column vector that includes a signal single
%   left side FFT;
% - a 2D matrix whose columns are signal single left
%   side FFTs;
% - a 3D matrix whose pages accommodate 2D matrix as
%    before.
% - f is the frequency vector (usually a column) where the
% FFTS have been evaluated.
% - fc is the center frequency for the triangFilter filter
% evaluation. fc can be both a scalar or a vector.
% - pcent is the percent of each fc central ferquencies.
% triangFilter smoothing function is a option to be used
% instead of KonnoOhmachiFilter.
%
% e.g. 
%      f = (0.1:0.1:10)'; x = f.^2; fc = (1:5)';  
%      y = triangFilter(x,f,fc,3);
%      semilogx(f,x,'b'); hold on; semilogx(fc,y,'r');
% e.g. 
%      f = (0.1:0.1:10)'; x = [f,f.^2]; fc = (1:5)';  
%      y = triangFilter(x,f,fc,3);
%      semilogx(f,x,'b'); hold on; semilogx(fc,y,'r'); 
% Written by D. Zuliani 2016/02/01

%
% MATRIX DIMENTION CHECK
nDims = ndims(x);
switch nDims
    case 2
        % NOTHING TO DO
    case 3
        % ACCOMODATING THE 3D STRUCTURE IN
        % A 2D MATRIX
        xSize = size(x);
        x = reshape(x,[xSize(1) xSize(2)*xSize(3) 1]);
    otherwise
end
%
% WORKING WITH COLUMN VECTORS
fc  = fc(:);
f   = f(:);
%
% GETTING CRUCIAL DIMs
m   = length(fc);   % number of KonnoHomachi center frequencies
n   = length(f);    % number of input frequencies
k   = size(x,2);    % number of input signals to filter
%
% Getting the signal percents
fInfBound=fc-fc*pcent/100;
fSupBound=fc+fc*pcent/100;

%
% Preallocating for speed up
%y       = zeros(1,length(fc));
%
% Filtering procedure
% for i=1:length(fc)
%     y(i)=mean(x(f>=fInfBound(i) & f<=fSupBound(i)).*triang(sum(f>=fInfBound(i) & f<=fSupBound(i))));
% end

%
% CREATING THE TRIANGULAR COEFFICIENTS.
% They are accommodated inside a NxM matrix: N=lenght(f), M=length(fc).
winTriang = zeros(length(f),length(fc)); % Preallocating for speed up
lenTriang = zeros(length(fc),1); % Preallocating for speed up
for i=1:length(fc)
    lenTriang(i) = sum(f>=fInfBound(i) & f<=fSupBound(i));
    winTriang((f>=fInfBound(i) & f<=fSupBound(i)),i)=triang(lenTriang(i));   
end
%
% CREATING A SINGLE COLUMN OF TRIANG COEFFS.
% All of them are clustered in groups fc1, fc2, ..., fcM. Each
% group is N length and includes  the traing coeffs. for
% a given central frequency fc
winTriang = reshape(winTriang,size(winTriang,2)*size(winTriang,1),1);
%
% REPEATING INPUT SIGNAL MATRIX ALONG M ROWS
% Each matrix is made by the splitted signals Sig1,Sig2,...,SigK
% vectors.
sig = repmat(x,m,1);
%
% ACCOMODATING THE COLUMNS OF x IN A SINGLE COLUMN
% The signal column vectors of the sig matrix are
% reshaped in a single column made by m*k*n rows.
% the 1st signal vector is made by the 1st n values,
% it is repeated m times and this is the 1st cluster
% including multiple copies of the 1st vector. The
% 2nd cluster is made in the same way for the 2nd
% signal.
sig = reshape(sig,m*k*n,1);
%
% REPEATING THE TRIANG COEFF. VECTOR.
% the winTriang vector is repeated to match the exact number
% of sig rows. The fc1,fc2,...,fcM groups are repeated
% k times along rows.
winTriang = repmat(winTriang,k,1);
%
% BUILDING THE SIMPLE sig*winTriang PRODUCT
% This is the product element by element between the
% signal vector and the column made by the fc groups.
% The operation is the 1st part of the triang
% filtering;
y = sig.*winTriang;
%
% WORKING WITH SUBINDEX
% They are needed to implement the full triang
% method using mean
subIdx = 1:1:m*k;
subIdx = repmat(subIdx,n,1);
subIdx = reshape(subIdx,n*m*k,1);
%
% TRIANG FILTERING LAST PART (mean values around the fc).
% The element by element products between the signal
% vector and the column made by the fc groups are finally
% summed group by group following the rule imposed by subIdx.
y = accumarray(subIdx,y)./repmat(lenTriang,k,1);
%
% RESHAPING
% y is a matrix made by k column vectors of m size.
y = (reshape(y,m,k));
%
% REBUILD A 3D MATRIX IF NEEDED
switch nDims
    case 2
        % NOTHING TO DO
    case 3
        % ACCOMODATING THE 2D STRUCTURE IN
        % A 3D MATRIX
        y = reshape(y,[m XSIZE(2) XSIZE(3)]);
    otherwise
end