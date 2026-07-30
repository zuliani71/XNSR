% Reproduce and validate the sanGiuliano_001 triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "sanGiuliano_001","Triang");
