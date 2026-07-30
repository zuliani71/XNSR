% Reproduce and validate the sanGiuliano_002 triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "sanGiuliano_002","Triang");
