function [xout] = ft2fft(xin,varargin);
switch length(varargin)
    case 0
        flag1 = 'o'; % it defines the KIND of simmetry, usually odd for real time series even/odd
        flag2 = 'o'; % it defines the KIND of the ending numbers of fft values (be careful) even/odd
    case 1
        flag1 = varargin{1};
    case 2
        flag1 = varargin{1};
        flag2 = varargin{2};
    otherwise
        ERMESSAGE='to many inputs';
        error(ERMESSAGE);
        xout=[];
        return;
end
xin=xin(:);
Points=size(xin,1);
if flag1=='e' %even
    xflip=transpose((fliplr(transpose(xin))));
elseif flag1=='o' %odd
    xflip=((fliplr(transpose(xin))))';
end
switch flag2
    case 'e'
        xout=[xin;xflip(2:end-1)];
    case 'o'
        xout=[xin;xflip(1:end-1)];
    otherwise
        ERMESSAGE='you must use a ''e'' or ''o'' flag as a 2nd parameter';
        error(ERMESSAGE);
        xout=[];
        return;
end