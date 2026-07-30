function [XN_DATA, outputPath] = generate_KonnoOhmachiCalibration(datasetName)
%GENERATE_KONNOOHMACHICALIBRATION Compatibility wrapper for old commands.
%   This function is retained for compatibility. New code should call:
%   update_XNSRCalibration(NAME,"KonnoOhmachi").

warning('XNSR:DeprecatedCalibrationGenerator', ...
    ['generate_KonnoOhmachiCalibration is retained for compatibility. ', ...
    'Use update_XNSRCalibration(NAME,"KonnoOhmachi") instead.']);
[XN_DATA, outputPath] = update_XNSRCalibration( ...
    string(datasetName),"KonnoOhmachi");
end
