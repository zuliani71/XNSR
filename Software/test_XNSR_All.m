function test_XNSR_All
%TEST_XNSR_ALL Run the complete portable XNSR regression suite.
%   After adding the cloned XNSR/Software directory to the MATLAB path,
%   call TEST_XNSR_ALL from any current working directory.

softwarePath = fileparts(mfilename('fullpath'));
testFiles = { ...
    'test_fft2ft.m', ...
    'test_readcfg.m', ...
    'test_triangFilter.m', ...
    'test_hv_konno.m', ...
    'test_KonnoOhmachiFilter.m', ...
    'test_XN_plotmatdata.m', ...
    'test_XN_Cruncher_Full.m'};

for testIndex = 1:numel(testFiles)
    fprintf('\n=== %s ===\n',testFiles{testIndex});
    runIsolated(fullfile(softwarePath,testFiles{testIndex}));
end
fprintf('\nTUTTA LA SUITE XNSR E'' SUPERATA\n');
end

function runIsolated(testPath)
% Each historical test is a script and may clear its own workspace.
run(testPath);
end
