% Made by D. Zuliani 2025/04/12
% This is a script to test XN_Cruncher on all the following datasets:
% 1) CA04
% 2) Edificio_Dorlando
% 3) Ferrara
% 4) Lorca_001
% 5) Lorca_002
% 6) Polignano;
% 7) sanGiuliano_001
% 8) sanGiuliano_002
% 9) ValMontanaia_001
%
% it will take different minutes depending on your CPU. You can try the
% test one by one if you want to reduce the time needed.
%
%% CLEANING
clear all;
close all;
format long g;
%
%% Running tests
test_XN_Cruncher_CA04;
test_XN_Cruncher_Edificio_Dorlando;
test_XN_Cruncher_Ferrara;
test_XN_Cruncher_Lorca_001;
test_XN_Cruncher_Lorca_002;
test_XN_Cruncher_Polignano;
test_XN_Cruncher_sanGiuliano_001;
test_XN_Cruncher_sanGiuliano_002;
test_XN_Cruncher_ValMontanaia_001;
close all;
clear all;