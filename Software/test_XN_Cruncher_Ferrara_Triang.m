% Reproduce and validate the Ferrara triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest("Ferrara","Triang");
