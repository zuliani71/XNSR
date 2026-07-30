% Reproduce and validate the Lorca_001 Konno-Ohmachi reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "Lorca_001","KonnoOhmachi");
