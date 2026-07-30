% Reproduce and validate the Lorca_002 triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest("Lorca_002","Triang");
