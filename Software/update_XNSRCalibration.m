function [XN_DATA, outputPath] = update_XNSRCalibration( ...
    datasetName,smoothingMethod)
%UPDATE_XNSRCALIBRATION Create or replace one official XNSR reference.
%   XN_DATA = UPDATE_XNSRCALIBRATION(NAME,METHOD) processes one of the nine
%   bundled datasets by calling XN_Cruncher and writes the result to the
%   corresponding DataCalib subdirectory.
%
%   METHOD can be "Triang" or "KonnoOhmachi". This maintenance function
%   overwrites an official reference and is not intended as a regression
%   test. Use test_XN_Cruncher_<dataset>_<method>.m to write to DataOut and
%   compare a new result without modifying DataCalib.

arguments
    datasetName (1,1) string
    smoothingMethod (1,1) string
end

[datasetName,fileList,repositoryPath] = ...
    get_XNSRDatasetFiles(datasetName);
smoothingMethod = validatestring(smoothingMethod, ...
    {'Triang','KonnoOhmachi'});

switch smoothingMethod
    case 'Triang'
        outputFamily = 'Triang';
        configName = 'XN_Cruncher_Triang.cfg';
        expectedSmoothingType = 'T';
    case 'KonnoOhmachi'
        outputFamily = 'KonnoOhmachi';
        configName = 'XN_Cruncher_Konno.cfg';
        expectedSmoothingType = 'K';
end

outputDirectory = fullfile(repositoryPath,'DataCalib',outputFamily);
configPath = fullfile(repositoryPath,'Cfg',configName);

assert(isfile(configPath), ...
    'XNSR:CalibrationUpdate:MissingConfig', ...
    'Calibration configuration file not found: %s',configPath);

if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
outputPath = fullfile(outputDirectory,[datasetName '.mat']);

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool;
end

calibrationTimer = tic;
XN_DATA = XN_Cruncher(fileList,outputPath,configPath,false);
elapsedSeconds = toc(calibrationTimer);

assert(isfile(outputPath), ...
    'XNSR:CalibrationUpdate:MissingOutput', ...
    'XN_Cruncher did not create %s.',outputPath);
assert(XN_DATA.PARAM.SCRIPT.SMOOTHING_WIN_TYPE == expectedSmoothingType, ...
    'XNSR:CalibrationUpdate:WrongSmoothingType', ...
    'XN_Cruncher used an unexpected smoothing method.');
if expectedSmoothingType == 'K'
    assert(XN_DATA.PARAM.KONNO.b == 40, ...
        'XNSR:CalibrationUpdate:WrongKonnoParameter', ...
        'The Konno-Ohmachi reference must use b = 40.');
else
    assert(XN_DATA.PARAM.TRIANG.perc == 3, ...
        'XNSR:CalibrationUpdate:WrongTriangParameter', ...
        'The triangular reference must use a 3%% window.');
end
assert(ismember(upper(XN_DATA.PARAM.SCRIPT.CALCULUS_MODE),{'M','MIXED'}));

numericFields = { ...
    'HV_RATIO','HV_RATIO_Fc','HV_STD','MAX_HV_RATIO','MAX_HV_F', ...
    'ALPHA_VEC','THETA_VEC','TXYZ'};
for fieldIndex = 1:numel(numericFields)
    fieldName = numericFields{fieldIndex};
    assert(isfield(XN_DATA,fieldName), ...
        'XNSR:CalibrationUpdate:MissingField', ...
        'XN_DATA.%s is missing.',fieldName);
    values = XN_DATA.(fieldName);
    assert(isnumeric(values) && ~isempty(values), ...
        'XNSR:CalibrationUpdate:InvalidField', ...
        'XN_DATA.%s must be a nonempty numeric array.',fieldName);
    assert(all(isfinite(values(:))), ...
        'XNSR:CalibrationUpdate:NonfiniteField', ...
        'XN_DATA.%s contains NaN or Inf.',fieldName);
end

alphaCount = floor(XN_DATA.PARAM.GEOM.MAX_ALPHA / ...
    XN_DATA.PARAM.GEOM.STEP_ALPHA) + 1;
thetaCount = floor(XN_DATA.PARAM.GEOM.MAX_THETA / ...
    XN_DATA.PARAM.GEOM.STEP_THETA) + 1;
orientationCount = alphaCount * thetaCount;
assert(size(XN_DATA.HV_RATIO,3) == orientationCount);
assert(numel(XN_DATA.MAX_HV_RATIO) == orientationCount);
assert(numel(XN_DATA.MAX_HV_F) == orientationCount);

fprintf(['RIFERIMENTO DATACALIB AGGIORNATO: %s | %s | ', ...
    '%d workers | %.3f s\n%s\n'], ...
    datasetName,outputFamily,pool.NumWorkers,elapsedSeconds,outputPath);
end
