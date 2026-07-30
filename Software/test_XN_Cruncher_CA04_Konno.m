% Full CA04 regression test using optimized Konno-Ohmachi smoothing.
% The output name is distinct from the historical CA04.mat calibration.

clearvars;
close all;
format long g;

softwarePath = fileparts(mfilename('fullpath'));
repositoryPath = fileparts(softwarePath);
dataPath = fullfile(repositoryPath,'DataIn');
outputPath = fullfile(repositoryPath,'DataOut', ...
    'CA04_Konno_optimized.mat');
configPath = fullfile(repositoryPath,'Cfg', ...
    'XN_Cruncher_Konno.cfg');

fileList = { ...
    fullfile(dataPath,'2004183110000.00.CA04.EHE.vel'), ...
    fullfile(dataPath,'2004183110000.00.CA04.EHN.vel'), ...
    fullfile(dataPath,'2004183110000.00.CA04.EHZ.vel')};

assert(all(cellfun(@isfile,fileList)),'CA04 input files are missing.');
assert(isfile(configPath),'The CA04 Konno configuration is missing.');

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool;
end
workerCount = pool.NumWorkers;

testTimer = tic;
XN_DATA = XN_Cruncher(fileList,outputPath,configPath);
totalSeconds = toc(testTimer);

assert(isfile(outputPath),'XN_Cruncher did not create the output MAT file.');
assert(isstruct(XN_DATA),'XN_Cruncher did not return a structure.');
calibrationPath = fullfile(repositoryPath,'DataCalib', ...
    'KonnoOhmachi','CA04.mat');
CALIBRATION_REPORT = compare_XNSRCalibration( ...
    outputPath,calibrationPath);
assert(XN_DATA.PARAM.SCRIPT.SMOOTHING_WIN_TYPE == 'K');
assert(ismember(upper(XN_DATA.PARAM.SCRIPT.CALCULUS_MODE),{'M','MIXED'}));

numericFields = { ...
    'HV_RATIO','HV_RATIO_Fc','HV_STD','MAX_HV_RATIO','MAX_HV_F', ...
    'ALPHA_VEC','THETA_VEC','TXYZ'};
for fieldIndex = 1:numel(numericFields)
    fieldName = numericFields{fieldIndex};
    assert(isfield(XN_DATA,fieldName), ...
        'XN_DATA is missing field %s.',fieldName);
    values = XN_DATA.(fieldName);
    assert(isnumeric(values) && ~isempty(values), ...
        'XN_DATA.%s must be a nonempty numeric array.',fieldName);
    assert(all(isfinite(values(:))), ...
        'XN_DATA.%s contains NaN or Inf.',fieldName);
end

alphaCount = floor(XN_DATA.PARAM.GEOM.MAX_ALPHA / ...
    XN_DATA.PARAM.GEOM.STEP_ALPHA) + 1;
thetaCount = floor(XN_DATA.PARAM.GEOM.MAX_THETA / ...
    XN_DATA.PARAM.GEOM.STEP_THETA) + 1;
orientationCount = alphaCount * thetaCount;

assert(numel(XN_DATA.ALPHA_VEC) == orientationCount);
assert(numel(XN_DATA.THETA_VEC) == orientationCount);
assert(size(XN_DATA.HV_RATIO,3) == orientationCount);
assert(numel(XN_DATA.MAX_HV_RATIO) == orientationCount);
assert(numel(XN_DATA.MAX_HV_F) == orientationCount);

fprintf('\nCA04 KONNO-OHMACHI COMPLETO SUPERATO\n');
fprintf('Worker paralleli:       %d\n',workerCount);
fprintf('Orientazioni elaborate: %d\n',orientationCount);
fprintf('Frequenze di uscita:    %d\n',numel(XN_DATA.HV_RATIO_Fc));
fprintf('Tempo complessivo:      %.3f s\n',totalSeconds);
fprintf('Output: %s\n',outputPath);
