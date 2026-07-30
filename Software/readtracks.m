function out=readtracks(varargin)
%READTRACKS Read a SAC or supported text waveform file.
%   OUT = READTRACKS(FILE) returns sampling frequency and waveform data.
%   OUT = READTRACKS(FILE1,FILE2,FILE3) reads three component files.
%
%   Originally written by D. Zuliani.
format long g;
% Initialize the output.
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
% Try SAC first, then fall back to supported text formats.
sacInfo=readsac(inputFile);
if ischar(sacInfo.SAC)
    % The file is likely an ASC, TXT, or TRC text waveform.
    textData = char(textread(inputFile,'%s','delimiter','\n'));
    i        = 1;
    iMax     = size(textData,1);
    stopFlag = 0;
    keyWords = {'Sampl. freq.:','Sampling rate:'};
    while (stopFlag==0 && i<iMax)
        if isempty(str2num(textData(i,:)))
            for j = 1:length(keyWords)
                if isempty(regexp(textData(i,:),keyWords{j}))
                else
                    out.samFreq = textscan(textData(i,length(keyWords{j})+1:end),'%f',1);
                    out.samFreq = out.samFreq{1};
                end
            end
        else
            dataStartIndex = i;
            stopFlag = 1;
        end
        i=i+1;
    end
    out.data=str2num(textData(dataStartIndex:iMax,:));
else
    % The file contains SAC data.
    out.samFreq = 1/sacInfo.Tsamp;
    out.data    = sacInfo.data(:,2);
end
