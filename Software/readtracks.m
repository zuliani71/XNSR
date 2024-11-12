function out=readtracks(varargin)
% Made by D. Zuliani 2016/02/02
format long g;
%
% DEFAULTS
out.samFreq = [];
out.data    = [];
switch nargin
    case 1
        inputFile = varargin{1};
    case 2
        disp('number of input arguments must be 1 or 3');
        return;
    case 3
        for i =1:nargin
            out{i}=readtracks(varargin{i});
        end
        return
end
%
% SAC or TXT
sacInfo=readsac(inputFile);
if ischar(sacInfo.SAC)
    % this is a txt file probably asc, txt  or trc type
    textData = char(textread(inputFile,'%s','delimiter','\n'));
    i        = 1;
    iMax     = size(textData,1);
    stopFlag = 0;
    keyWords = {'Sampl. freq.:'};
    while (stopFlag==0 && i<iMax)
        if isempty(str2num(textData(i,:)))
            if isempty(regexp(textData(i,:),keyWords{1}))
            else
                out.samFreq = textscan(textData(i,length(keyWords{1})+1:end),'%f',1);
                out.samFreq = out.samFreq{1};
            end
        else
            dataStartIndex = i;
            stopFlag = 1;
        end
        i=i+1;
    end
    out.data=str2num(textData(dataStartIndex:iMax,:));
else
    % this is a sac file
    out.samFreq = 1/sacInfo.Tsamp;
    out.data    = sacInfo.data(:,2);    
end
