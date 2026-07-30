% Reproduce and validate the Lorca_001 triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest("Lorca_001","Triang");
