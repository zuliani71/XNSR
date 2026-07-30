% Reproduce and validate the CA04 triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest("CA04","Triang");
