
function HV = hv_specgram(XYZ,varargin)
% HV = hv_specgram(XYZ,varagin)calculates
% the H/V ratio over the XYZ dataset.
% The calculus is performed using the specgram
% matlab function.
% XYZ is a matrix of column vectors.
% 1st vector is the N component
% 2nd vector is the E component
% 3rd vector is the V component
% By default the sampling rate is 1s = 100Hz, you can change it
%   using a 3rd parameter
% By default the moving windows size used by hv_specgram is 50, you can change
%   it using a 4th parameter
% By default the fft size used by hv_specgram is 256, you can change it using
%   a 5th parameter
% By default the output is given at 1Hz step frequencies from 0 to the half of the
%   sampling rate, you can change it using a 6th parameter
%
% e.g. hv_specgram(XYZ,100,50,256,(0:0.01:50)) will work at 100Hz with 50
% samples of moving window overlaps, an fft size of 256 and an output
% frequency vector defined by 0 to 50Hz with 0.01Hz steps

%
% DEFAULTS
SAMPFREQ=   100;
WINSIZE =   50;
FFTSIZE =   256;
FREQLIST=   (1:1:SAMPFREQ/2)';
%
% DEALING WITH INPUT ARGUMENTS
switch length(varargin)
    case 1
        if ~isempty(varargin{1})
            SAMPFREQ=varargin{1};
        end
    case 2
        if ~isempty(varargin{1})
            SAMPFREQ=varargin{1};
        end
        if ~isempty(varargin{2})
            WINSIZE=varargin{2};
        end        
    case 3
        if ~isempty(varargin{1})
            SAMPFREQ=varargin{1};
        end
        if ~isempty(varargin{2})
            WINSIZE=varargin{2};
        end
        if ~isempty(varargin{3})
            FFTSIZE=varargin{3};
        end
    case 4
        if ~isempty(varargin{1})
            SAMPFREQ=varargin{1};
        end
        if ~isempty(varargin{2})
            WINSIZE=varargin{2};
        end
        if ~isempty(varargin{3})
            FFTSIZE=varargin{3};
        end
        if ~isempty(varargin{4})
            FREQLIST=varargin{4};
        end
end
%
% DEALING WITH INPUT DATASET
H1  =   XYZ(:,1);
H2  =   XYZ(:,2);
V   =   XYZ(:,3);
%
% DOING GEOMETRIC AVERAGE
H   =   sqrt(H1.*H2);
%
% WORKING ON THE VERTICAL COMPONENT
[SV,~,~] = specgram(V,FREQLIST,SAMPFREQ,FFTSIZE,WINSIZE);
SV=mean(abs(SV),2);
%
% WORKING ON THE HORIZONTAL COMPONENTS
[SH,FH, ~] = specgram(H,FREQLIST,SAMPFREQ,FFTSIZE,WINSIZE);
SH=mean(abs(SH),2);
%
% H/V ratio
HV = [FH,SH./SV];