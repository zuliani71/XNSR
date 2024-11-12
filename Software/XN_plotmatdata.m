function XN_DATA=XN_plotmatdata(varargin)
% Made by D. Zuliani 2013/09/19
%

%
% DEFAULTS
%close all;
format long g;
scrsz           =   get(0,'ScreenSize');
FONT.SIZE       =   16;
FONT.WEIGHT     =   'Bold';
FONT.NAME       =   'Courier';
XN_DATA         =   [];
MAINPLOTTYPE    =   '3D';
%
% PARSING INPUT ARGUMENTS
switch length(varargin)
    case 1
        DEFPATH = varargin{1};
    otherwise
        DEFPATH = '/Users/dzuliani/SHARED/Projects/Seismology/2013.HV_RATIO_XNSR/DATI/DAROSARIAGALLIPOLI/Archivio';
end
%
% Load Matlab XN dataset
[FILENAME,PATHNAME] = uigetfile('*.mat','Select the MATLAB XN Dataset',DEFPATH);
XN_DATA=load([PATHNAME,FILENAME]);
XN_DATA=XN_DATA.XN_DATA;
%
% PLOTS
%
% PRELIMINAR SIGNAL PLOTS
figure('Position',[1 1 scrsz(3)*0.365 scrsz(4)/3])
%
% TIME DOMAIN PLOTS
subplot(3,2,1);
H(1) = plot(XN_DATA.TXYZ(:,1),XN_DATA.TXYZ(:,2),'r');
title('TIME DOMAIN');
xlabel('T(s)');
ylabel('X(V)');
grid on;
axis tight
subplot(3,2,3);
H(3) = plot(XN_DATA.TXYZ(:,1),XN_DATA.TXYZ(:,3),'b');
xlabel('T(s)');
ylabel('Y(V)');
grid on;
axis tight
subplot(3,2,5);
H(5) = plot(XN_DATA.TXYZ(:,1),XN_DATA.TXYZ(:,4),'m');
xlabel('T(s)');
ylabel('Z(V)');
grid on;
axis tight
%
% FREQUENCY DOMAIN PLOTS
XFFT = fft2ft(abs(fft(XN_DATA.TXYZ(:,2))),XN_DATA.PARAM.SIG.F);
YFFT = fft2ft(abs(fft(XN_DATA.TXYZ(:,3))),XN_DATA.PARAM.SIG.F);
ZFFT = fft2ft(abs(fft(XN_DATA.TXYZ(:,4))),XN_DATA.PARAM.SIG.F);
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
% FILTERED X/N RATIO PLOTS
figure('Position',[1 scrsz(4)/2 scrsz(3)*0.365 scrsz(4)/2])
XN_DATA.SUB_PLT2=subplot(1,2,2);
plot(1:10,1:10);
set(gca,...
    'FontSize',FONT.SIZE,...
    'FontWeight',FONT.WEIGHT,...
    'FontName',FONT.NAME);
MINF    = min(XN_DATA.MAX_HV_F);
MAXF    = max(XN_DATA.MAX_HV_F);
STEPF   = (MAXF-MINF)/10;
ICOLOR  = XN_DATA.MAX_HV_F;
IDIM    = exp(6*(XN_DATA.MAX_HV_RATIO/max(XN_DATA.MAX_HV_RATIO))); % max H/V modulus proportional to circle radius
I       = find(IDIM==0);
IDIM(I)=1;
axes('position', [0.05,0.1,0.4,0.8]);
switch MAINPLOTTYPE
    case {'2D','2'}
        % HQ=scatter(180/pi*XN_DATA.ALPHA_VEC,180/pi*XN_DATA.THETA_VEC,IDIM,ICOLOR,'filled'); %old matlab
        HQ=scatter(180/pi*XN_DATA.ALPHA_VEC,180/pi*XN_DATA.THETA_VEC,IDIM(:),ICOLOR(:),'filled'); %working on matlab R2018a for mac
        %set(gca,'YDir','reverse');
        axis ij;
        %axis tight;
        xlabel('AZIMUTH ANGLE [degrees]',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        ylabel('DIP ANGLE [degrees]',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        axis([-5,XN_DATA.PARAM.GEOM.MAX_ALPHA+5,-5,XN_DATA.PARAM.GEOM.MAX_THETA+5]);
        grid on;
    case {'3D','3'}
        H=stem3(180/pi*XN_DATA.ALPHA_VEC,180/pi*XN_DATA.THETA_VEC,XN_DATA.MAX_HV_RATIO,'color','k');
        set(H,'Marker','none');
        hold on;
        % scatter3(180/pi*XN_DATA.ALPHA_VEC,180/pi*XN_DATA.THETA_VEC,XN_DATA.MAX_HV_RATIO,IDIM(:),ICOLOR(:),'filled'); %old matlab
        scatter3(180/pi*XN_DATA.ALPHA_VEC,180/pi*XN_DATA.THETA_VEC,XN_DATA.MAX_HV_RATIO(:,:),IDIM(:),ICOLOR(:),'filled'); %working on matlab R2018a for mac
        xlabel('AZIMUTH ANGLE [degrees]',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        ylabel('DIP ANGLE [degrees]',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        zlabel('X/N RATIO',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        axis([-5,XN_DATA.PARAM.GEOM.MAX_ALPHA+5,-5,XN_DATA.PARAM.GEOM.MAX_THETA+5,0,10*mean(XN_DATA.MAX_HV_RATIO)]);
        grid on;
end
title('MAX(X/N) with X/N=func(AZIMUTH,DIP,f)',...
    'FontSize',FONT.SIZE,...
    'FontWeight',FONT.WEIGHT,...
    'FontName',FONT.NAME);
%
% COLORBAR X FREQUENCY
YTICK_VEC = 0:STEPF:MAXF;
YTICK_LAB = (cellstr(num2str(YTICK_VEC')))';
YTICK_VEC = (YTICK_VEC/MAXF)*256;
H=colorbar('FontSize',FONT.SIZE,...
    'FontWeight',FONT.WEIGHT,...
    'FontName',FONT.NAME);
title(H,'f[Hz]',...
    'FontSize',FONT.SIZE,...
    'FontWeight',FONT.WEIGHT,...
    'FontName',FONT.NAME);
set(gca,...
    'FontSize',FONT.SIZE,...
    'FontWeight',FONT.WEIGHT,...
    'FontName',FONT.NAME);
%
% SAVING THE MAIN DATASET INFOS INSIDE THE FIGURE HANDLE
set(gca,'UserData',XN_DATA);
%
% BUILDING THE CURSOR MODE FEATURE
POINTEROBJ = datacursormode;
set (POINTEROBJ,'Enable','on',...
    'DisplayStyle','datatip',...
    'UpdateFcn',@doratio);
%
% FUNCTION doratio FOR "ON THE FLY" SPECTRAL RATIO CALCULUS
    function output_txt = doratio(obj,event_obj)
        % Display the position of the data cursor
        % obj          Currently not used (empty)
        % event_obj    Handle to event object
        % RATIO        Works with the H/V ratios.
        XN_DATA = get(gca,'UserData');
        pos = get(event_obj,'Position');
%         output_txt = {['AZIMUTH ANGLE: ',num2str(pos(1),4),'°'],...
%             ['DIP ANGLE: ',num2str(pos(2),4),'°']};
        output_txt = {['AZIMUTH=',num2str(pos(1),4),'°'],...
            ['DIP=',num2str(pos(2),4),'°']};
        %
        % RECOVER SELECTED ALPHA and THETA
        I1 = find(round(XN_DATA.ALPHA_VEC/pi*180)==round(pos(1)));
        I2 = find(round(XN_DATA.THETA_VEC/pi*180)==round(pos(2)));
        I = intersect(I1,I2);
        THETA   =   XN_DATA.THETA_VEC(I)/pi*180;
        ALPHA   =   XN_DATA.ALPHA_VEC(I)/pi*180;
        %
        % WORKING WITH DATA CALCULATED BY MATRIX MANIPULATIONS
%         output_txt{end+1} = ['X/N MAX RATIO:           ',num2str(XN_DATA.MAX_HV_RATIO(I),4)];
%         output_txt{end+1} = ['X/N MAX RATIO freq:(Hz): ',num2str(XN_DATA.MAX_HV_F(I),4)];
%         output_txt{end+1} = ['MAX(X/N)=',num2str(XN_DATA.MAX_HV_RATIO(I),4)];
%         output_txt{end+1} = ['f_MAX(X/N)=',num2str(XN_DATA.MAX_HV_F(I),4),'Hz'];
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
        set(gca,...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        xlabel('f(Hz)',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        ylabel('X/N ratio',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        TITLE_STRING = ['AZIMUTH=',num2str(pos(1),4),'° ',...
            'DIP=',num2str(pos(2),4),'° ',...
            'MAX(X/N)=',num2str(XN_DATA.MAX_HV_RATIO(I),4),' ',...
            'f_{MAX(X/N)}=',num2str(XN_DATA.MAX_HV_F(I),4),'Hz'];
        title(TITLE_STRING,...
            'Interpreter','tex',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        
        %         axis ([min(XN_DATA.HV_RATIO_Fc),...
        %             max(XN_DATA.HV_RATIO_Fc),...
        %             0,20]);
        axis ([min(XN_DATA.HV_RATIO_Fc),...
            max(XN_DATA.HV_RATIO_Fc),...
            min(min(XN_DATA.HV_RATIO-XN_DATA.HV_STD)),...
            max(max(XN_DATA.HV_RATIO+XN_DATA.HV_STD))]);
        legend([PLT_HV,PLT_STDP],{'MEAN X/N RATIO','STD DEVIATION'},...
            'Location','NorthWest',...
            'FontSize',FONT.SIZE,...
            'FontWeight',FONT.WEIGHT,...
            'FontName',FONT.NAME);
        hold off;
    end
end