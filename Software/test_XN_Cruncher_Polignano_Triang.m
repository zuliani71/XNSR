% Reproduce and validate the Polignano triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest("Polignano","Triang");
