function [XN_DATA,CALIBRATION_REPORT,outputPath] = ...
    run_XNSRTest(datasetName,smoothingMethod)
%RUN_XNSRTEST Run one bundled experiment and validate it against DataCalib.
%   RUN_XNSRTEST(NAME,"Triang") writes DataOut/Triang/NAME.mat.
%   RUN_XNSRTEST(NAME,"KonnoOhmachi") writes
%   DataOut/KonnoOhmachi/NAME.mat.

arguments
    datasetName (1,1) string
    smoothingMethod (1,1) string
end

% Make the test reproducible even when another project provides functions
% with the same legacy names (for example readsac, readtracks or fft2ft).
softwarePath = fileparts(mfilename('fullpath'));
addpath(softwarePath,'-begin');

[datasetName,fileList,repositoryPath] = ...
    get_XNSRDatasetFiles(datasetName);
smoothingMethod = string(validatestring(smoothingMethod, ...
    {'Triang','KonnoOhmachi'}));

switch smoothingMethod
    case "Triang"
        outputFamily = 'Triang';
        configName = 'XN_Cruncher_Triang.cfg';
    case "KonnoOhmachi"
        outputFamily = 'KonnoOhmachi';
        configName = 'XN_Cruncher_Konno.cfg';
end

configPath = fullfile(repositoryPath,'Cfg',configName);
outputDirectory = fullfile(repositoryPath,'DataOut',outputFamily);
calibrationPath = fullfile(repositoryPath,'DataCalib', ...
    outputFamily,datasetName + ".mat");

assert(isfile(configPath), ...
    'XNSR:Test:MissingConfig', ...
    'Test configuration file not found: %s',configPath);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
outputPath = fullfile(outputDirectory,datasetName + ".mat");

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool;
end

testTimer = tic;
XN_DATA = XN_Cruncher(fileList,outputPath,configPath,false);
elapsedSeconds = toc(testTimer);
CALIBRATION_REPORT = compare_XNSRCalibration( ...
    outputPath,calibrationPath);

fprintf(['TEST XNSR SUPERATO: %s | %s | %d workers | %.3f s\n', ...
    'Output: %s\n'],datasetName,outputFamily,pool.NumWorkers, ...
    elapsedSeconds,outputPath);
end
