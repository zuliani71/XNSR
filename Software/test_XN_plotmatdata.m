% Regression tests for XN_plotmatdata.

softwarePath = fileparts(mfilename('fullpath'));
repositoryPath = fileparts(softwarePath);
calibrationPath = fullfile(repositoryPath,'DataCalib','Triang','CA04.mat');
assert(isfile(calibrationPath));

originalVisibility = get(groot,'DefaultFigureVisible');
cleanup = onCleanup(@() restoreFigures(originalVisibility));
set(groot,'DefaultFigureVisible','off');
close all force;

data2D = XN_plotmatdata(calibrationPath,'2D');
assert(isstruct(data2D));
assert(numel(findall(groot,'Type','figure')) == 2);
close all force;

data3D = XN_plotmatdata(calibrationPath,'3D');
assert(isstruct(data3D));
assert(numel(findall(groot,'Type','figure')) == 2);
close all force;

temporaryDirectory = tempname;
mkdir(temporaryDirectory);
temporaryCleanup = onCleanup(@() rmdir(temporaryDirectory,'s'));

loaded = load(calibrationPath,'XN_DATA');
XN_DATA = loaded.XN_DATA;
XN_DATA.MAX_HV_RATIO(1) = NaN;
XN_DATA.MAX_HV_F(1) = NaN;
nanPath = fullfile(temporaryDirectory,'nan_dataset.mat');
save(nanPath,'XN_DATA');
nanData = XN_plotmatdata(nanPath,'2D');
assert(isstruct(nanData));
close all force;

XN_DATA.MAX_HV_RATIO(:) = NaN;
allNanPath = fullfile(temporaryDirectory,'all_nan_dataset.mat');
save(allNanPath,'XN_DATA');
assertError(@() XN_plotmatdata(allNanPath), ...
    'XNSR:PlotMatData:NoFiniteOrientations');
assert(isempty(findall(groot,'Type','figure')));

invalidPath = fullfile(temporaryDirectory,'invalid.mat');
invalidContent = 1;
save(invalidPath,'invalidContent');
assertError(@() XN_plotmatdata(invalidPath), ...
    'XNSR:PlotMatData:InvalidDataset');
assertError(@() XN_plotmatdata(calibrationPath,'4D'), ...
    'XNSR:PlotMatData:InvalidPlotType');

clear temporaryCleanup cleanup
fprintf('TUTTI I TEST XN_PLOTMATDATA SONO SUPERATI\n');

function assertError(callback,expectedIdentifier)
failedAsExpected = false;
try
    callback();
catch exception
    failedAsExpected = strcmp(exception.identifier,expectedIdentifier);
end
assert(failedAsExpected, ...
    'Expected error %s was not raised.',expectedIdentifier);
end

function restoreFigures(originalVisibility)
close all force;
set(groot,'DefaultFigureVisible',originalVisibility);
end
