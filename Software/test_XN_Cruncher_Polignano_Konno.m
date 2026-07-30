% Reproduce and validate the Polignano Konno-Ohmachi reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "Polignano","KonnoOhmachi");
