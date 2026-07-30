% Reproduce and validate all nine triangular references.
clearvars;
close all;
format long g;
datasetNames = [ ...
    "CA04","Edificio_Dorando","Ferrara","Lorca_001","Lorca_002", ...
    "Polignano","sanGiuliano_001","sanGiuliano_002", ...
    "ValMontanaia_001"];
for datasetIndex = 1:numel(datasetNames)
    run_XNSRTest(datasetNames(datasetIndex),"Triang");
    close all;
end
