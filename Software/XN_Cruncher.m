function XN_DATA=XN_Cruncher(varargin)
% Made by D. Zuliani and T. Stabile 2013/07/18
% Modified by D. Zuliani 2013/07/19
% Modified by D. Zuliani 2013/08/20
% Modified by D. Zuliani 2013/08/22
% Modified by D. Zuliani 2013/08/23
% Modified by D. Zuliani 2013/08/24
% Modified by D. Zuliani 2013/08/27
% Modified by D. Zuliani 2013/09/19
% Modified by D. Zuliani 2013/09/20
% Modified by D. Zuliani 2016/02/02

%
% 1st TIME REMEBER TO OPEN THE MATLABPOOL
% matlabpool open
%
% CLEANING
clear all;
close all;
format long g;
%
% DEFAULTS PARAMETERS
scrsz = get(0,'ScreenSize');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PARAMETERS YOU CAN CHANGE STARTS HERE %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GEMETRIC PARAMETERS
PARAM.GEOM.MAX_ALPHA       = 180;          %MAX AZIMUTH ANGLE, MIN = 0 BY DEFAULT
PARAM.GEOM.MAX_THETA       = 30;           %MAX DIP ANGLE, MIN = 0 BY DEFAULT
PARAM.GEOM.STEP_ALPHA      = 10;           %AZIMUTH ANGLE DEGREE STEPS
PARAM.GEOM.STEP_THETA      = 2;            %DIP ANGLE DEGREE STEPS
%
% SIGNAL PARAMETERS
%PARAM.SIG.F               = 128;          %DATALOGGING FREQUENCY
PARAM.SIG.FFTSIZE         = 3840;         %COMMON FFT SIZE (confortable to use for all methods)
PARAM.SIG.T_LIM           = [0,600];    %TIME INTERVAL OF INPUT SIGNAL TO USE (TEST 600s and 900s)
PARAM.SIG.F_LIM           = [0.2,10];     %OUTPUT SIGNAL FREQUENCY INTERVAL TO PLOT
%
% PREFILTERING (BUTTEERWORTH) AND TAPERING (TUKEY WINDOW) PARAMS
TUKEYWINPARAM   = 0.05;         %Tapering ratio applied to main signal
Fn              = [0.5,20];     %band bass filter frequency limits e.g. [0.1,25];
Or              = 4;            %band pass filter order
%
% Common Method PARAMETERS:
PARAM.COMM.fstp= 0.01; %OUTPUT FREQUENCY STEP (Hz)
PARAM.COMM.fc  = (PARAM.SIG.F_LIM(1):PARAM.COMM.fstp:PARAM.SIG.F_LIM(2))'; %OUTPUT FREQUENCY VECTOR
%
% KonnoOhmachi PARAMETERS:
PARAM.KONNO.b   = 40; % b KonnoOhmachi parameter
%
% Triang PARAMETERS:
PARAM.TRIANG.perc= 3; % +/-3% of signala at every fc frequency
%
% SCRIPT CONTROL PARAMETERS
DOTAPERING          = 'N';
DOFILTERING         = 'Y';
DODETRENDING        = 'Y';
MAINPLOTTYPE        = '2D'; %Available values are 3D and 2D
SMOOTHING_WIN_TYPE  = 'T';  %T for Triangulare windows, K for Konnomachi method
CALCULUS_MODE       = 'M';  %Available methods are 'V' for Full Vectorization, or 'M' for partial loop and vectorization
%
% INPUT DATA
% DATAPATH        = '/Users/dzuliani/VBOX.SHARE/Projects/Seismology/2013.HV_RATIO/DATI/DAROSARIAGALLIPOLI/Archivio';
% FILELIST        = {'Edificio_Dorando.asc'};
% DATAPATH        = '/Users/dzuliani/VBOX.SHARE/Projects/Seismology/2013.HV_RATIO/DATI/2013.CNR.TONY_STABILE/Archivio';
% FILELIST        = {'EW.002','NS.001','Z.003'};
DATAPATH = '/Users/dzuliani/SHARED/Projects/Seismology/2013.HV_RATIO_XNSR/DATI/CRS';
FILELIST = {'2004183110000.00.CA04.EHE.vel','2004183110000.00.CA04.EHN.vel','2004183110000.00.CA04.EHZ.vel'};
if length(FILELIST) == 1
    CURRENTFILE = [DATAPATH,'/',FILELIST{1}];
    dataSig=readtracks(CURRENTFILE);
    if size(dataSig.data,2) < 3
        disp('Input file does not include all the components needed');
        return
    else
        Y           = dataSig.data(:,1);
        X           = dataSig.data(:,2);
        Z           = dataSig.data(:,3);
        PARAM.SIG.F = dataSig.samFreq;  %DATALOGGING FREQUENCY
    end
else
    for i = 1:length(FILELIST)
        CURRENTFILE = [DATAPATH,'/',FILELIST{i}];
        dataSig=readtracks(CURRENTFILE);
        switch i
            case 1
                Y           = dataSig.data;
                PARAM.SIG.F = dataSig.samFreq; %DATALOGGING FREQUENCY
            case 2
                X           = dataSig.data;
                PARAM.SIG.F = dataSig.samFreq; %DATALOGGING FREQUENCY
            case 3
                Z           = dataSig.data;
                PARAM.SIG.F = dataSig.samFreq; %DATALOGGING FREQUENCY
            otherwise
                disp ('FILE NOT FOUND')
                return
        end
    end
end
%
% OUTPUT DATA
FILE_MATLAB_OUT =   [DATAPATH,'/','XN_Ratio.mat'];


%
% TIME VECTOR
T = (0:1/PARAM.SIG.F:(1/PARAM.SIG.F)*(length(X)-1));
I = find((T<=PARAM.SIG.T_LIM(2)) & (T>=PARAM.SIG.T_LIM(1)));
T = T(I)';
X = X(I)';
Y = Y(I)';
Z = Z(I)';
%
% PRELIMINAR FILTERING (BUTTERWORTH)
switch DOFILTERING
    case {'Y','y','Yes','yes','YES'}
        Wcs     =   2*pi*PARAM.SIG.F;   %Carrier beat frequency
        Wcn     =   2*pi*Fn/(Wcs/2);    %Normalized cutoff frequency.
        [B,A]   =   butter(Or,Wcn);     %Band Pass Butterworth Filter
        X       =   filtfilt(B,A,X);    %Filtered signal (zero-phase digital filtering)
        Y       =   filtfilt(B,A,Y);    %Filtered signal (zero-phase digital filtering)
        Z       =   filtfilt(B,A,Z);    %Filtered signal (zero-phase digital filtering)
    otherwise
end
%
% SIGNAL TAPERING
switch DOTAPERING
    case {'Y','y','Yes','yes','YES'}
        TAPESIZE = length(X);
        X=X.*(tukeywin(TAPESIZE,TUKEYWINPARAM))';
        Y=Y.*(tukeywin(TAPESIZE,TUKEYWINPARAM))';
        Z=Z.*(tukeywin(TAPESIZE,TUKEYWINPARAM))';
    otherwise
end
%
% SIGNAL DETRENDING
switch DODETRENDING
    case {'Y','y','Yes','yes','YES'}
        X = detrend(X,'constant');
        Y = detrend(Y,'constant');
        Z = detrend(Z,'constant');
    otherwise
end
%
% PRELIMINAR SIGNAL PLOTS
figure('Position',[1 1 scrsz(3)*0.365 scrsz(4)/3])
%
% TIME DOMAIN PLOTS
subplot(3,2,1);
H(1) = plot(T,X,'r');
title('TIME DOMAIN');
xlabel('T(s)');
ylabel('X(V)');
grid on;
axis tight
subplot(3,2,3);
H(3) = plot(T,Y,'b');
xlabel('T(s)');
ylabel('Y(V)');
grid on;
axis tight
subplot(3,2,5);
H(5) = plot(T,Z,'m');
xlabel('T(s)');
ylabel('Z(V)');
grid on;
axis tight
%
% FREQUENCY DOMAIN PLOTS
XFFT = fft2ft(abs(fft(X)),PARAM.SIG.F);
YFFT = fft2ft(abs(fft(Y)),PARAM.SIG.F);
ZFFT = fft2ft(abs(fft(Z)),PARAM.SIG.F);
subplot(3,2,2);
H(2) = semilogx(XFFT(:,1),20*log10(XFFT(:,2)),'r');
title('FREQUENCY DOMAIN');
xlabel('f(Hz)');
ylabel('X(dB)');
grid on;
axis tight
subplot(3,2,4);
H(4) = semilogx(YFFT(:,1),20*log10(YFFT(:,2)),'b');
xlabel('f(Hz)');
ylabel('Y(dB)');
grid on;
axis tight
subplot(3,2,6);
H(6) = semilogx(ZFFT(:,1),20*log10(ZFFT(:,2)),'m');
xlabel('f(Hz)');
ylabel('Z(dB)');
grid on;
axis tight
%
% SIGNAL SPLIT
disp('SIGNAL SPLIT');
tic;
X_SPLIT = detrend(WinSplit(X,PARAM.SIG.FFTSIZE,PARAM.SIG.FFTSIZE*0.1));
Y_SPLIT = detrend(WinSplit(Y,PARAM.SIG.FFTSIZE,PARAM.SIG.FFTSIZE*0.1));
Z_SPLIT = detrend(WinSplit(Z,PARAM.SIG.FFTSIZE,PARAM.SIG.FFTSIZE*0.1));
toc
%
% 1) WORKING INSIDE THE FFT DOMAIN
% 2) RECOVERING THE FFT HALF LEFT SIDE
disp('FFT + FFT2FT');
tic;
FT_X    = fft2ft(fft(X_SPLIT),PARAM.SIG.F);
FT_Y    = fft2ft(fft(Y_SPLIT),PARAM.SIG.F);
FT_Z    = fft2ft(fft(Z_SPLIT),PARAM.SIG.F);
toc
%
% REDUCING THE DATASET ACCORDING THE FREQUENCY LIMITS
disp('DASET REDUCING BY FREQ. LIMS');
tic;
F_VECT      = fft2ft(fft(X_SPLIT(:,1)),PARAM.SIG.F);
F_VECT      = F_VECT(:,1);
I           = find((F_VECT>=PARAM.SIG.F_LIM(1)) & (F_VECT<=PARAM.SIG.F_LIM(2)));
F_VECT      = F_VECT(I);
FT_X        = FT_X(I,:,2);
FT_Y        = FT_Y(I,:,2);
FT_Z        = FT_Z(I,:,2);
toc
%
%
%%%%%%%%%%%%%%% MANIPULATION MATRIX STARTS HERE %%%%%%%%%%%%%%%
disp('DATASET ROTATIONS');
tic;
%
% AZIMUTH ANGLE
ALPHA = 0:PARAM.GEOM.STEP_ALPHA:PARAM.GEOM.MAX_ALPHA;
ALPHA = ALPHA/180*pi;
%
% DIP ANGLE
THETA = 0:PARAM.GEOM.STEP_THETA:PARAM.GEOM.MAX_THETA;
THETA = THETA/180*pi;
%
% BUILDING ANGLES MATRIX
ALPHA_VEC   = repmat(ALPHA',length(THETA),1);
THETA_VEC   = reshape(repmat(THETA',1,length(ALPHA))',length(ALPHA_VEC),1);
%
% preallocation for speeding up
FT_X_ROT   = zeros([size(FT_X),size(ALPHA_VEC,1)]);
FT_Y_ROT   = FT_X_ROT;
FT_Z_ROT   = FT_X_ROT;
for i=1:size(ALPHA_VEC,1)
    FT_X_ROT(:,:,i)   = cos(THETA_VEC(i))*cos(ALPHA_VEC(i))*FT_X - cos(THETA_VEC(i))*sin(ALPHA_VEC(i))*FT_Y + sin(ALPHA_VEC(i))*FT_Z;
    FT_Y_ROT(:,:,i)   = sin(ALPHA_VEC(i))*FT_X + cos(ALPHA_VEC(i))*FT_Y;
    FT_Z_ROT(:,:,i)   = -sin(THETA_VEC(i))*cos(ALPHA_VEC(i))*FT_X + sin(ALPHA_VEC(i))*sin(THETA_VEC(i))*FT_Y + cos(THETA_VEC(i))*FT_Z;
end
toc
%
% HORIZONTAL COMPONENTS MEAN
disp('HORIZONTAL COMPONENTS MEAN');
tic;
%FT_XY_ROT   = abs(FT_X_ROT+FT_Y_ROT)/2;
FT_XY_ROT   = sqrt(abs(FT_X_ROT).*abs(FT_Y_ROT));
FT_Z_ROT    = abs(FT_Z_ROT);
toc
% XN_DATA = FT_XY_ROT;
% return
%
switch upper(SMOOTHING_WIN_TYPE)
    case 'K'
        % KONNOOMACHI FILTERING
        disp('WORKING WITH KONNOOHMACHI');
        switch upper(CALCULUS_MODE)
            case {'MIXED','M'}
                %
                % preallocation for speeding up
                XY_SPECTRUM =   zeros(size(PARAM.COMM.fc,1),size(FT_XY_ROT,2),size(FT_XY_ROT,3));
                Z_SPECTRUM  =   XY_SPECTRUM;
                %
                % SLOWER FOR A SMALL AMMOUNT OF MEMORY USAGE BUT
                % FASTER WITH BIG DASASETS
                tic
                size(FT_XY_ROT,3)
                parfor i = 1:size(FT_XY_ROT,3)
                    XY_SPECTRUM(:,:,i) = KonnoOhmachiFilter(FT_XY_ROT(:,:,i),F_VECT,PARAM.COMM.fc,PARAM.KONNO.b);
                    Z_SPECTRUM(:,:,i)  = KonnoOhmachiFilter(FT_Z_ROT(:,:,i),F_VECT,PARAM.COMM.fc,PARAM.KONNO.b);
                end
                toc
            case {'VECTORIZATION','V'}
                %
                % SLOWER FOR BIG AMMOUNT OF MEMORY USAGE BUT
                % FASTER WITH SMALL DASASETS
                XY_SPECTRUM = KonnoOhmachiFilter(FT_XY_ROT(:,:,:),F_VECT,PARAM.COMM.fc,PARAM.KONNO.b);
                Z_SPECTRUM   = KonnoOhmachiFilter(FT_Z_ROT(:,:,:),F_VECT,PARAM.COMM.fc,PARAM.KONNO.b);
                XY_SPECTRUM = squeeze(XY_SPECTRUM);
                Z_SPECTRUM  = squeeze(Z_SPECTRUM);
        end
    case 'T'
        % TRIANGULAR FILTERING
        disp('WORKING WITH TRIANGULAR SMOOTHING');
        switch upper(CALCULUS_MODE)
            case {'MIXED','M'}
                %
                % preallocation for speeding up
                XY_SPECTRUM =   zeros(size(PARAM.COMM.fc,1),size(FT_XY_ROT,2),size(FT_XY_ROT,3));
                Z_SPECTRUM  =   XY_SPECTRUM;
                %
                % SLOWER FOR A SMALL AMMOUNT OF MEMORY USAGE BUT
                % FASTER WITH BIG DASASETS
                tic
                parfor i = 1:size(FT_XY_ROT,3)
                    XY_SPECTRUM(:,:,i) = triangFilter(FT_XY_ROT(:,:,i),F_VECT,PARAM.COMM.fc,PARAM.TRIANG.perc);
                    Z_SPECTRUM(:,:,i)  = triangFilter(FT_Z_ROT(:,:,i),F_VECT,PARAM.COMM.fc,PARAM.TRIANG.perc);
                end
                toc
            case {'VECTORIZATION','V'}
                %
                % SLOWER FOR BIG AMMOUNT OF MEMORY USAGE BUT
                % FASTER WITH SMALL DASASETS
                XY_SPECTRUM = triangFilter(FT_XY_ROT(:,:,:),F_VECT,PARAM.COMM.fc,PARAM.TRIANG.perc);
                Z_SPECTRUM  = triangFilter(FT_Z_ROT(:,:,:),F_VECT,PARAM.COMM.fc,PARAM.TRIANG.perc);
                XY_SPECTRUM = squeeze(XY_SPECTRUM);
                Z_SPECTRUM  = squeeze(Z_SPECTRUM);
        end
        
end
%
% H/V ratio
HV_RATIO                = XY_SPECTRUM./Z_SPECTRUM;
HV_STD                  = std(HV_RATIO,1,2);
HV_RATIO                = mean(HV_RATIO,2);
%HV_STD                  = (1./Z_SPECTRUM).*sqrt(XY_STD.^2+(XY_SPECTRUM.^2)./(Z_SPECTRUM.^2).*(Z_STD.^2));
%
% RECOVERING MAX MODULE FREQUENCY RESPONSES ABD THEIR FREQUENCY VALUES
[MAX_HV_RATIO,I]        = max(HV_RATIO,[],1);
MAX_HV_F                = PARAM.COMM.fc(I);
XN_DATA=[];
%%%%%%%%%%%%%%% MANIPULATION MATRIX STOPS HERE %%%%%%%%%%%%%%%
%
% PLOTS
figure('Position',[1 scrsz(4)/2 scrsz(3)*0.365 scrsz(4)/2])
SUB_PLT2=subplot(1,2,2);
plot(1:10,1:10);
MINF    = min(MAX_HV_F);
MAXF    = max(MAX_HV_F);
STEPF   = (MAXF-MINF)/10;
ICOLOR  = MAX_HV_F;
IDIM    = exp(6*(MAX_HV_RATIO/max(MAX_HV_RATIO))); % max H/V modulus proportional to circle radius
I       = find(IDIM==0);
IDIM(I)=1;
axes('position', [0.05,0.1,0.4,0.8]);
switch MAINPLOTTYPE
    case {'2D','2'}
        scatter(180/pi*ALPHA_VEC,180/pi*THETA_VEC,IDIM(:),ICOLOR(:),'filled');
        axis ij;
        xlabel('AZIMUTH ANGLE [degrees]')
        ylabel('DIP ANGLE [degrees]')
        axis([-5,PARAM.GEOM.MAX_ALPHA+5,-5,PARAM.GEOM.MAX_THETA+5]);
        grid on;
    case {'3D','3'}
        H=stem3(180/pi*ALPHA_VEC,180/pi*THETA_VEC,MAX_HV_RATIO,'color','k');
        set(H,'Marker','none');
        hold on;
        scatter3(180/pi*ALPHA_VEC,180/pi*THETA_VEC,MAX_HV_RATIO(:),IDIM(:),ICOLOR','filled');
        xlabel('AZIMUTH ANGLE [degrees]')
        ylabel('DIP ANGLE [degrees]')
        zlabel('X/N RATIO');
        axis([-5,PARAM.GEOM.MAX_ALPHA+5,-5,PARAM.GEOM.MAX_THETA+5,0,10*mean(MAX_HV_RATIO)]);
        grid on;
end
title('X/N Method');
%
% COLORBAR X FREQUENCY
YTICK_VEC = 0:STEPF:MAXF;
YTICK_LAB = (cellstr(num2str(YTICK_VEC')))';
YTICK_VEC = (YTICK_VEC/MAXF)*256;
H=colorbar;
title(H,'f[Hz]');
%
% WORKING WITH SELECTIONS
XN_DATA.HV_RATIO    =   HV_RATIO;
XN_DATA.HV_RATIO_Fc =   PARAM.COMM.fc;
XN_DATA.HV_STD      =   HV_STD;
XN_DATA.MAX_HV_RATIO=   MAX_HV_RATIO;
XN_DATA.MAX_HV_F    =   MAX_HV_F;
XN_DATA.ALPHA_VEC   =   ALPHA_VEC;
XN_DATA.THETA_VEC   =   THETA_VEC;
XN_DATA.SUB_PLT2    =   SUB_PLT2;
XN_DATA.TXYZ        =   [T,X',Y',Z'];
XN_DATA.PARAM       =   PARAM;
%
% SAVE THE DATASET
save(FILE_MATLAB_OUT,'XN_DATA');
%
% SAVING THE MAIN DATASET INFOS INSIDE THE FIGURE HANDLE
set(gca,'UserData',XN_DATA);
%
% BUILDING THE CURSOR MODE FEATURE
POINTEROBJ = datacursormode;
set (POINTEROBJ,'Enable','on',...
    'DisplayStyle','window',...
    'UpdateFcn',@doratio);
%
% FUNCTIONS
%
% FUNCTION doratio FOR "ON THE FLY" SPECTRAL RATIO CALCULUS
    function output_txt = doratio(obj,event_obj)
        % Display the position of the data cursor
        % obj          Currently not used (empty)
        % event_obj    Handle to event object
        % RATIO        Works with the H/V ratios.
        XN_DATA = get(gca,'UserData');
        pos = get(event_obj,'Position');
        output_txt = {['AZIMUTH ANGLE: ',num2str(pos(1),4),'°'],...
            ['DIP ANGLE: ',num2str(pos(2),4),'°']};
        %
        % RECOVER SELECTED ALPHA and THETA
        I1 = find(round(XN_DATA.ALPHA_VEC/pi*180)==round(pos(1)));
        I2 = find(round(XN_DATA.THETA_VEC/pi*180)==round(pos(2)));
        I = intersect(I1,I2);
        THETA   =   XN_DATA.THETA_VEC(I)/pi*180;
        ALPHA   =   XN_DATA.ALPHA_VEC(I)/pi*180;
        %
        % WORKING WITH DATA CALCULATED BY MATRIX MANIPULATIONS
        output_txt{end+1} = ['X/N MAX RATIO:           ',num2str(XN_DATA.MAX_HV_RATIO(I),4)];
        output_txt{end+1} = ['X/N MAX RATIO freq:(Hz): ',num2str(XN_DATA.MAX_HV_F(I),4)];
        subplot(XN_DATA.SUB_PLT2);
        PLT_HV=semilogx(XN_DATA.HV_RATIO_Fc,...
            XN_DATA.HV_RATIO(:,I));
        hold on;
        PLT_STDP=semilogx(XN_DATA.HV_RATIO_Fc,...
            XN_DATA.HV_RATIO(:,I)+XN_DATA.HV_STD(:,I),'r');
        PLT_STDM=semilogx(XN_DATA.HV_RATIO_Fc,...
            XN_DATA.HV_RATIO(:,I)-XN_DATA.HV_STD(:,I),'r');
        semilogx(XN_DATA.MAX_HV_F(I),...
            XN_DATA.MAX_HV_RATIO(I),...
            'bo','MarkerFaceColor','b');
        grid on;
        xlabel('f(Hz)');
        ylabel('X/N ratio');
        axis ([min(XN_DATA.HV_RATIO_Fc),...
            max(XN_DATA.HV_RATIO_Fc),...
            min(min(XN_DATA.HV_RATIO-XN_DATA.HV_STD)),...
            max(max(XN_DATA.HV_RATIO+XN_DATA.HV_STD))]);
        legend([PLT_HV,PLT_STDP],{'MEAN X/N RATIO','STD DEVIATION'})
        hold off;
    end
end