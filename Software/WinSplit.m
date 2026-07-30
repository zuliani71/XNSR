function Y = WinSplit(X,varargin)
%WINSPLIT Split signals into overlapping, tapered windows.
%   Y = WINSPLIT(X) uses an eight-sample Hamming window with 50% overlap.
%   A vector X produces a matrix whose columns are windowed segments.
%   A matrix X is interpreted as one signal per row and produces a
%   three-dimensional array with one page per input signal.
%
%   Y = WINSPLIT(X,WIN) uses WIN directly when it is a vector. When WIN is
%   an integer, it specifies the segment length and a Hamming window of
%   that length is generated.
%
%   Y = WINSPLIT(X,WIN,N) uses an overlap of N samples. N must be smaller
%   than the window length. The default overlap is 50%.
%
%   DIMS = WINSPLIT(X,WIN,N,'D') returns the output dimensions without
%   computing the segments, which is useful for preallocation.
%
%   Originally written by D. Zuliani.

% Defaults.
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
% Resolve the requested window.
switch length(WIN)
    case 1
        WIN = bartlett(WIN);
    otherwise
end
WIN=WIN(:);
% Normalize input orientation to one signal per column.
if size(X,1)==1 && size(X,2)>1
    % Convert row vectors to columns.
    X = X(:);    
elseif size(X,2)==1 && size(X,1)>1
    % Column vectors already have the required orientation.
else
    % Matrix rows contain signals; transpose them into columns.
    X=X.';
end
% Compute the number and shape of output segments.
WIN_SIZE= length(WIN);
NROWS   = WIN_SIZE;     % Number of output rows.
SIGDIM  = size(X,1);    % Samples per input signal.
NSIGS   = size(X,2);    % Number of input signals.
NCOLS   = fix((SIGDIM-WIN_OVERLAP)/(WIN_SIZE-WIN_OVERLAP)); % Output segments.
% Build segment indices.
switch upper(WIN_MODE)
    case {'D','DIM','DIMS','DIMENTION','DIMENTIONS'}
        Y       = [SIGDIM,NSIGS,NCOLS];
    case {'S','STD','STANDARD'}
        ICOLS   = 1+(0:(NCOLS-1))*(WIN_SIZE-WIN_OVERLAP);
        IROWS   = (1:NROWS)';
        if size(X,2) > 1
            % Extract each segment.
            Y = X(IROWS(:,ones(1,NCOLS))+ICOLS(ones(WIN_SIZE,1),:)-1,:);
            % Apply the selected taper.
            WIN=(reshape(repmat(WIN,1,NCOLS),1,WIN_SIZE*NCOLS))';
            Y= reshape(bsxfun(@times,Y,WIN),WIN_SIZE,NCOLS,NSIGS);        
        else
            % Extract each segment.
            Y = X(IROWS(:,ones(1,NCOLS))+ICOLS(ones(WIN_SIZE,1),:)-1);
            % Apply the selected taper.
            Y= bsxfun(@times,Y,WIN);
        end
end
