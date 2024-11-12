
function HV = hv_simplefilter(XYZ,varargin)
% HV = hv_simplefilter(XYZ,varagin)calculates
% the H/V ratio over the XYZ dataset.
% The calculus is performed using the
% a simple fft ratio and a smoothing
% function of a moving averaging window.
% XYZ is a matrix of column vectors.
% 1st vector is the N component
% 2nd vector is the E component
% 3rd vector is the V component
% By default the sampling rate is 1s = 100Hz, you can change it
%   using a 3rd parameter
% By default the moving window size parameter is empty (no filtering),
%   you can change it using a 4th parameter
% By default the fft size used by hv_simplefilter is 256, you can
%   change it using a 5th parameter
%
% e.g. hv_simplefilter(XYZ,100,50,512); will work at 100Hz os sampling
% frequency with a moving window average of size 50 samples and a fft
% size of 512

% DEFAULTS
SAMPFREQ=   100;
WINSIZE  =   [];
FFTSIZE =   256;
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
SV = fft2ft(abs(fft(V,FFTSIZE)),SAMPFREQ);
% WORKING ON THE HORIZONTAL COMPONENTS
SH = fft2ft(abs(fft(H,FFTSIZE)),SAMPFREQ);
%
% WORKING WITH THE FREQUENCY VECTOR
FH = SH(:,1);
%
% MOVING WINDOW AVERAGE SMOOTHING
if isempty(WINSIZE)
    SH = SH(:,2:end);
    SV = SV(:,2:end);
else
    SH = filter(ones(1,WINSIZE)/WINSIZE,1,SH(:,2:end));
    SV = filter(ones(1,WINSIZE)/WINSIZE,1,SV(:,2:end));
end
%
% H/V ratio
HV = [FH,SH./SV];