function XN_DATA = XN_plotmatdata(varargin)
%XN_PLOTMATDATA Display an XNSR MAT dataset.
%   XN_PLOTMATDATA() selects a MAT file from DataOut using a 2D view.
%   XN_PLOTMATDATA(PATH) accepts either a MAT file or an initial folder.
%   XN_PLOTMATDATA(VIEW) selects the default folder and VIEW ('2D'/'3D').
%   XN_PLOTMATDATA(PATH,VIEW) specifies both source and view.
%
%   Cancelling the file-selection dialog returns an empty array.

format long g;
font.Size = 16;
font.Weight = 'Bold';
font.Name = 'Courier';
XN_DATA = [];
plotType = '2D';

softwarePath = fileparts(mfilename('fullpath'));
dataSource = fullfile(softwarePath,'..','DataOut');
narginchk(0,2);
if nargin == 1
    candidate = upper(char(string(varargin{1})));
    if ismember(candidate,{'2D','2','3D','3'})
        plotType = candidate;
    else
        dataSource = char(string(varargin{1}));
    end
elseif nargin == 2
    dataSource = char(string(varargin{1}));
    plotType = upper(char(string(varargin{2})));
end

assert(ismember(plotType,{'2D','2','3D','3'}), ...
    'XNSR:PlotMatData:InvalidPlotType', ...
    'Plot type must be ''2D'' or ''3D''.');
assert(isfile(dataSource) || isfolder(dataSource), ...
    'XNSR:PlotMatData:MissingSource', ...
    'MAT file or initial folder not found: %s',dataSource);

if isfile(dataSource)
    matPath = dataSource;
else
    [fileName,pathName] = uigetfile( ...
        '*.mat','Select the MATLAB XN Dataset',dataSource);
    if isequal(fileName,0) || isequal(pathName,0)
        return
    end
    matPath = fullfile(pathName,fileName);
end

matVariables = whos('-file',matPath);
assert(any(strcmp({matVariables.name},'XN_DATA')), ...
    'XNSR:PlotMatData:InvalidDataset', ...
    '%s does not contain an XN_DATA variable.',matPath);
loaded = load(matPath,'XN_DATA');
assert(isstruct(loaded.XN_DATA) && ...
    isscalar(loaded.XN_DATA), ...
    'XNSR:PlotMatData:InvalidDataset', ...
    '%s does not contain a scalar XN_DATA structure.',matPath);
XN_DATA = loaded.XN_DATA;
validateDataset(XN_DATA,matPath);

frequency = XN_DATA.HV_RATIO_Fc(:);
alphaDegrees = 180/pi*XN_DATA.ALPHA_VEC(:);
thetaDegrees = 180/pi*XN_DATA.THETA_VEC(:);
maximumRatio = XN_DATA.MAX_HV_RATIO(:);
maximumFrequency = XN_DATA.MAX_HV_F(:);
orientationCount = numel(alphaDegrees);
ratio = reshape(XN_DATA.HV_RATIO,numel(frequency),orientationCount);
ratioStd = reshape(XN_DATA.HV_STD,numel(frequency),orientationCount);

validOrientation = isfinite(alphaDegrees) & isfinite(thetaDegrees) & ...
    isfinite(maximumRatio) & isfinite(maximumFrequency);
assert(any(validOrientation), ...
    'XNSR:PlotMatData:NoFiniteOrientations', ...
    'XN_DATA contains no finite orientation results to plot.');

screenSize = get(groot,'ScreenSize');
figure('Position',[1 1 screenSize(3)*0.365 screenSize(4)/3]);
plotSignalAndSpectrum(XN_DATA);

figure('Position', ...
    [1 screenSize(4)/2 screenSize(3)*0.365 screenSize(4)/2]);
ratioAxes = subplot(1,2,2);
plot(ratioAxes,NaN,NaN);
set(ratioAxes,'FontSize',font.Size, ...
    'FontWeight',font.Weight,'FontName',font.Name);

plotAlpha = alphaDegrees(validOrientation);
plotTheta = thetaDegrees(validOrientation);
plotRatio = maximumRatio(validOrientation);
plotFrequency = maximumFrequency(validOrientation);
positiveMaximum = max(plotRatio);
if positiveMaximum > 0
    markerSize = exp(6*max(plotRatio,0)/positiveMaximum);
else
    markerSize = 36*ones(size(plotRatio));
end

mainAxes = axes('Position',[0.05,0.1,0.4,0.8]);
switch plotType
    case {'2D','2'}
        scatter(mainAxes,plotAlpha,plotTheta,markerSize, ...
            plotFrequency,'filled');
        set(mainAxes,'YDir','reverse');
    case {'3D','3'}
        stem3(mainAxes,plotAlpha,plotTheta,plotRatio, ...
            'Color','k','Marker','none');
        hold(mainAxes,'on');
        scatter3(mainAxes,plotAlpha,plotTheta,plotRatio, ...
            markerSize,plotFrequency,'filled');
        zlabel(mainAxes,'X/N RATIO', ...
            'FontSize',font.Size,'FontWeight',font.Weight, ...
            'FontName',font.Name);
end
xlabel(mainAxes,'AZIMUTH ANGLE [degrees]', ...
    'FontSize',font.Size,'FontWeight',font.Weight,'FontName',font.Name);
ylabel(mainAxes,'DIP ANGLE [degrees]', ...
    'FontSize',font.Size,'FontWeight',font.Weight,'FontName',font.Name);
title(mainAxes,'MAX(X/N) with X/N=func(AZIMUTH,DIP,f)', ...
    'FontSize',font.Size,'FontWeight',font.Weight,'FontName',font.Name);
grid(mainAxes,'on');
axis(mainAxes,'tight');
colorHandle = colorbar(mainAxes);
title(colorHandle,'f[Hz]', ...
    'FontSize',font.Size,'FontWeight',font.Weight,'FontName',font.Name);
set(mainAxes,'FontSize',font.Size, ...
    'FontWeight',font.Weight,'FontName',font.Name);

XN_DATA.SUB_PLT2 = ratioAxes;
set(mainAxes,'UserData',XN_DATA);
pointerObject = datacursormode(ancestor(mainAxes,'figure'));
set(pointerObject,'Enable','on','DisplayStyle','datatip', ...
    'UpdateFcn',@displayRatio);

    function outputText = displayRatio(~,eventObject)
        position = get(eventObject,'Position');
        [~,orientationIndex] = min( ...
            (alphaDegrees-position(1)).^2 + ...
            (thetaDegrees-position(2)).^2);
        outputText = { ...
            ['AZIMUTH=',num2str(alphaDegrees(orientationIndex),4),'°'], ...
            ['DIP=',num2str(thetaDegrees(orientationIndex),4),'°']};

        meanRatio = ratio(:,orientationIndex);
        standardDeviation = ratioStd(:,orientationIndex);
        validFrequency = isfinite(frequency) & frequency > 0 & ...
            isfinite(meanRatio) & isfinite(standardDeviation);
        axes(ratioAxes);
        cla(ratioAxes);
        if ~any(validFrequency)
            text(ratioAxes,0.5,0.5,'No finite spectral ratio available', ...
                'HorizontalAlignment','center');
            return
        end

        ratioLine = semilogx(ratioAxes,frequency(validFrequency), ...
            meanRatio(validFrequency));
        hold(ratioAxes,'on');
        deviationLine = semilogx(ratioAxes,frequency(validFrequency), ...
            meanRatio(validFrequency)+standardDeviation(validFrequency), ...
            'r');
        semilogx(ratioAxes,frequency(validFrequency), ...
            meanRatio(validFrequency)-standardDeviation(validFrequency), ...
            'r');
        if isfinite(maximumFrequency(orientationIndex)) && ...
                isfinite(maximumRatio(orientationIndex))
            semilogx(ratioAxes,maximumFrequency(orientationIndex), ...
                maximumRatio(orientationIndex), ...
                'bo','MarkerFaceColor','b');
        end
        grid(ratioAxes,'on');
        axis(ratioAxes,'tight');
        xlabel(ratioAxes,'f(Hz)');
        ylabel(ratioAxes,'X/N ratio');
        title(ratioAxes,sprintf( ...
            ['AZIMUTH=%.4g° DIP=%.4g° MAX(X/N)=%.4g ', ...
             'f_{MAX(X/N)}=%.4g Hz'], ...
            alphaDegrees(orientationIndex), ...
            thetaDegrees(orientationIndex), ...
            maximumRatio(orientationIndex), ...
            maximumFrequency(orientationIndex)), ...
            'Interpreter','tex');
        legend(ratioAxes,[ratioLine,deviationLine], ...
            {'MEAN X/N RATIO','STD DEVIATION'}, ...
            'Location','NorthWest');
        set(ratioAxes,'FontSize',font.Size, ...
            'FontWeight',font.Weight,'FontName',font.Name);
        hold(ratioAxes,'off');
    end
end

function validateDataset(data,matPath)
requiredFields = {'TXYZ','PARAM','HV_RATIO','HV_RATIO_Fc','HV_STD', ...
    'MAX_HV_RATIO','MAX_HV_F','ALPHA_VEC','THETA_VEC'};
missing = requiredFields(~isfield(data,requiredFields));
assert(isempty(missing), ...
    'XNSR:PlotMatData:MissingFields', ...
    'XN_DATA in %s is missing: %s',matPath,strjoin(missing,', '));

validateattributes(data.TXYZ,{'numeric'}, ...
    {'2d','nonempty','ncols',4},mfilename,'XN_DATA.TXYZ');
assert(isfield(data.PARAM,'SIG') && isfield(data.PARAM.SIG,'F') && ...
    isnumeric(data.PARAM.SIG.F) && isscalar(data.PARAM.SIG.F) && ...
    isfinite(data.PARAM.SIG.F) && data.PARAM.SIG.F > 0, ...
    'XNSR:PlotMatData:InvalidSamplingFrequency', ...
    'XN_DATA.PARAM.SIG.F must be a positive finite scalar.');

numericFields = requiredFields(~strcmp(requiredFields,'PARAM'));
for fieldIndex = 1:numel(numericFields)
    fieldName = numericFields{fieldIndex};
    assert(isnumeric(data.(fieldName)) && ~isempty(data.(fieldName)), ...
        'XNSR:PlotMatData:InvalidField', ...
        'XN_DATA.%s must be a nonempty numeric array.',fieldName);
end

orientationCount = numel(data.ALPHA_VEC);
assert(numel(data.THETA_VEC) == orientationCount && ...
    numel(data.MAX_HV_RATIO) == orientationCount && ...
    numel(data.MAX_HV_F) == orientationCount, ...
    'XNSR:PlotMatData:OrientationSizeMismatch', ...
    'Orientation vectors and maximum-result vectors have different sizes.');
frequencyCount = numel(data.HV_RATIO_Fc);
assert(numel(data.HV_RATIO) == frequencyCount*orientationCount && ...
    numel(data.HV_STD) == frequencyCount*orientationCount, ...
    'XNSR:PlotMatData:RatioSizeMismatch', ...
    'HV_RATIO and HV_STD dimensions do not match frequency/orientation data.');
end

function plotSignalAndSpectrum(data)
time = data.TXYZ(:,1);
colors = {'r','b','m'};
labels = {'X','Y','Z'};
for component = 1:3
    subplot(3,2,2*component-1);
    plot(time,data.TXYZ(:,component+1),colors{component});
    if component == 1
        title('TIME DOMAIN');
    end
    xlabel('T(s)');
    ylabel([labels{component},'(V)']);
    grid on;
    axis tight

    spectrum = fft2ft(abs(fft(data.TXYZ(:,component+1))), ...
        data.PARAM.SIG.F);
    subplot(3,2,2*component);
    semilogx(spectrum(:,1), ...
        20*log10(max(spectrum(:,2),realmin)),colors{component});
    if component == 1
        title('FREQUENCY DOMAIN');
    end
    xlabel('f(Hz)');
    ylabel([labels{component},'(dB)']);
    grid on;
    axis tight
end
end
