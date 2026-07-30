% Reproduce and validate the CA04 Konno-Ohmachi reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest("CA04","KonnoOhmachi");
