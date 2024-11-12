function Y = WinSplit(X,varargin)
% 
% Y = WinSplit(X,varargin) divides X into eight
% segments with 50% overlap, each segment is
% windowed with a Hamming window.
% X can be both a signal vector (column or row) or
% a signal matrix with signals distributed by rows.
%
% If X is a single vector (column or row), Y is a
% bidimentional matrix and each segment produced by
% WinSplit is accommodated inside the columns of Y.
% If X is a matrix (made of signals distributed by
% rows), Y is a 3D matrix organzied by pages. Each
% page includes a 2D matrix made by columns of
% WinSplit segments. WinSplit yelds one page for
% every X row.
%
% Y = WinSplit(X,WIN) when WIN is a vector, divides
% X into segments of length equal to the length of
% WIN, and then windows each segment with the vector
% specified by WIN.
% If WIN is an integer, X is divided into segments
% of length equal to that integer value, and a
% Hamming window of equal length is used.  If win 
% is not specified, the default (Hamming window of
% length 8) is used.
% Y = WinSplit(X,WIN,N) N is the number of samples
% each segment of X overlaps. N must be an integer
% smaller than WIN if WIN is an integer. N must be
% an integer smaller than the length of WIN if WIN
% is a vector. If N is not specified, the default
% value is used to obtain a 50% overlap.
%
% DIMS = WinSplit(X,WIN,N,MODE) gives out the dimentions
% of Y avoiding calculus when MODE ='D'. That is
% confortable when you want to preallocate memory before
% running WinSplit(X,WIN).
%
% Made by D. Zuliani 2013/08/19

%
% DEFAULTS
WIN         = round(length(X)/8);
WIN_OVERLAP = round(WIN*0.5);
WIN_MODE    = 'STD';
switch length(varargin)
    case 1
        if ~isempty(varargin{1})
            WIN         = varargin{1};
            WIN_OVERLAP = round(WIN*0.5);
        end
    case 2
        if ~isempty(varargin{1})
            WIN         = varargin{1};
        end
        if ~isempty(varargin{2})
            WIN_OVERLAP = varargin{2};
        else
            WIN_OVERLAP = round(WIN*0.5);
        end
    case 3
        if ~isempty(varargin{1})
            WIN         = varargin{1};
        end
        if ~isempty(varargin{2})
            WIN_OVERLAP = varargin{2};
        else
            WIN_OVERLAP = round(WIN*0.5);
        end
        if ~isempty(varargin{3})
            WIN_MODE = varargin{3};
        end
end
%
%DEALING WITH WIN
switch length(WIN)
    case 1
        WIN = bartlett(WIN);
    otherwise
end
WIN=WIN(:);
%
% DEALING WITH INPUT VECTOR
if size(X,1)==1 && size(X,2)>1
    % is X a row? Change it to a column
    X = X(:);    
elseif size(X,2)==1 && size(X,1)>1
    % is X a column? Leave it as it is
else
    % X is a matrix with signals organized by rows
    % so a change to a matrix with signals organized
    % by columns is neeeded
    X=X.';
end
%
% DEALING WITH OUTPUT VECTOR
WIN_SIZE= length(WIN);
NROWS   = WIN_SIZE;     %number of output rows
SIGDIM  = size(X,1);    %size of signal included in the X matrix
NSIGS   = size(X,2);    %number of signals included inside the X matrix
NCOLS   = fix((SIGDIM-WIN_OVERLAP)/(WIN_SIZE-WIN_OVERLAP));    %number of output segments
%
% INDEXING
switch upper(WIN_MODE)
    case {'D','DIM','DIMS','DIMENTION','DIMENTIONS'}
        Y       = [SIGDIM,NSIGS,NCOLS];
    case {'S','STD','STANDARD'}
        ICOLS   = 1+(0:(NCOLS-1))*(WIN_SIZE-WIN_OVERLAP);
        IROWS   = (1:NROWS)';
        if size(X,2) > 1
            %
            % BUILD Y UP
            Y = X(IROWS(:,ones(1,NCOLS))+ICOLS(ones(WIN_SIZE,1),:)-1,:);
            %
            % TAPPERING
            WIN=(reshape(repmat(WIN,1,NCOLS),1,WIN_SIZE*NCOLS))';
            Y= reshape(bsxfun(@times,Y,WIN),WIN_SIZE,NCOLS,NSIGS);        
        else
            %
            % BUILD Y UP
            Y = X(IROWS(:,ones(1,NCOLS))+ICOLS(ones(WIN_SIZE,1),:)-1);
            %
            % TAPPERING
            Y= bsxfun(@times,Y,WIN);
        end
end