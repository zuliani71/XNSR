% Reproduce and validate the sanGiuliano_001 Konno-Ohmachi reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "sanGiuliano_001","KonnoOhmachi");
