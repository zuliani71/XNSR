
function HV = hv_konno(XYZ,varargin)
% HV = hv_konno(XYZ,varagin)calculates
% the H/V ratio over the XYZ dataset.
% The calculus is performed using the
% a simple fft ratio and a smoothing
% function by KonnoOhmachi.
% XYZ is a matrix of column vectors.
% 1st vector is the N component
% 2nd vector is the E component
% 3rd vector is the V component
% By default the sampling rate is 1s = 100Hz, you can change it
%   using a 3rd parameter
% By default the KonnoOhmachi paremeter used by hv_konno is 40, you can change
%   it using a 4th parameter
% By default the fft size used by hv_konno is the full length of the signal,
%   you can change it using a 5th parameter
%
%
% e.g. hv_konno(XYZ,100,40,1024); will work at 100Hz os sampling frequency
% with the KonnoOhmachi parameter set at 40 and the fft size set to 1024

% DEFAULTS
SAMPFREQ =   100;
KONNOPAR =   40;
FFTSIZE  =   512;
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
            KONNOPAR=varargin{2};
        end        
    case 3
        if ~isempty(varargin{1})
            SAMPFREQ=varargin{1};
        end
        if ~isempty(varargin{2})
            KONNOPAR=varargin{2};
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
SV = fft(V,FFTSIZE)/FFTSIZE;
SV = 2*abs(SV(1:ceil(FFTSIZE/2)));
% WORKING ON THE HORIZONTAL COMPONENTS
SH = fft(H,FFTSIZE)/FFTSIZE;
SH = 2*abs(SH(1:ceil(FFTSIZE/2)));
%
% WORKING WITH THE FREQUENCY VECTOR
FH = (SAMPFREQ/2*linspace(0,1,ceil(FFTSIZE/2)))';
% 
% KonnoOhmachi SMOOTHING
SV  = KonnoOhmachi(SV,FH,KONNOPAR);
SH  = KonnoOhmachi(SH,FH,KONNOPAR);
%
% H/V ratio
HV = [FH,(SH./SV)];