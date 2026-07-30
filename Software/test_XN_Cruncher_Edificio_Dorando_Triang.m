% Reproduce and validate the Edificio_Dorando triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "Edificio_Dorando","Triang");
