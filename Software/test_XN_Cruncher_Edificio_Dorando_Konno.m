% Reproduce and validate the Edificio_Dorando Konno-Ohmachi reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "Edificio_Dorando","KonnoOhmachi");
