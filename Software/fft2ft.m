function [xout] = fft2ft(xin,varargin);
% xout = fft2ft(xin) implements the fft to ft (single left side fft)
% conversion between the xin and xout vectors. If 2nd parameter is
% inserted this last is intended to be the sample rate frequency.
% In this case xout is a vector with 2 collumns 1st for frequency
% values 2nd for amplitude values.
% fft2ft now can handle 2D and 3D input matrix. For each matrix
% signal vectors are accommdated inside the matrix collumns.
% If 2nd parameter is inserted this last is intended to be the sample
% rate frequency. In this case xout is a 3D matrix with at least 2
% pages, the 1st for frequency values and the 2nd and further for
% amplitude values.
% e.g.
% xout = fft2ft(xin); % no sampling rate defined
% xout = fft2ft(xout,500); %500Hz sampling rate
%
% Made by D. Zuliani 2004
% Modified by D. Zuliani 2013/07/19
% Modified by D. Zuliani 2013/08/19

%
% CHECK OUT IF INPUT IS A VECTOR OR A MATRIX
DIM = size(xin);
if isempty(find(DIM==1,1))
    Points=size(xin,1);
    switch ndims(xin)
        case 2
            % 2D MATRIX
            % WARNING THIS WORKS JUST ALONG COLLUMNS
            if mod(Points,2)==0
                xout=xin(1:1+Points/2,:);
            else
                xout=xin(1:(1+Points)/2,:);
            end
        case 3
            % 3D MATRIX
            % WARNING THIS WORKS JUST ALONG COLLUMNS AND PAGES
            if mod(Points,2)==0
                xout=xin(1:1+Points/2,:,:);
            else
                xout=xin(1:(1+Points)/2,:,:);
            end
    end
    %
    % WORKING WITH FREQUENCY MATRIX
    if nargin == 2
        fsample = varargin{1};
        f = ((0:1:(size(xout,1)-1))*(fsample/2)/size(xout(:,1),1))';
        f = repmat(f,1,size(xout,2));
        xout=(2/Points)*xout;
        xout=cat(3,f,xout);
    else
    end
else
    %
    % VECTOR
    xin=xin(:);
    Points=size(xin,1);
    if mod(Points,2)==0
        xout=xin(1:1+Points/2);
    else
        xout=xin(1:(1+Points)/2);
    end
    if nargin == 2
        fsample = varargin{1};
        f = ((0:1:(size(xout,1)-1))*(fsample/2)/size(xout(:,1),1))';
        xout=2/Points*xout;
        xout=[f,xout];
    else
    end
end