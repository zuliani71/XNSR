% Reproduce and validate all nine experiments with both smoothing methods.
clearvars;
close all;
format long g;
datasetNames = [ ...
    "CA04","Edificio_Dorando","Ferrara","Lorca_001","Lorca_002", ...
    "Polignano","sanGiuliano_001","sanGiuliano_002", ...
    "ValMontanaia_001"];
smoothingMethods = ["Triang","KonnoOhmachi"];
for methodIndex = 1:numel(smoothingMethods)
    for datasetIndex = 1:numel(datasetNames)
        run_XNSRTest(datasetNames(datasetIndex), ...
            smoothingMethods(methodIndex));
        close all;
    end
end
