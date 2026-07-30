% Reproduce and validate the Ferrara Konno-Ohmachi reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "Ferrara","KonnoOhmachi");
