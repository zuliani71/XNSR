% Reproduce and validate the ValMontanaia_001 triangular reference.
clearvars;
close all;
format long g;
[XN_DATA,CALIBRATION_REPORT] = run_XNSRTest( ...
    "ValMontanaia_001","Triang");
