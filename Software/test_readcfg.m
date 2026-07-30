% Regression tests for the XNSR configuration parser.

softwarePath = fileparts(mfilename('fullpath'));
repositoryPath = fileparts(softwarePath);

triang = readcfg(fullfile(repositoryPath,'Cfg', ...
    'XN_Cruncher_Triang.cfg'));
konno = readcfg(fullfile(repositoryPath,'Cfg', ...
    'XN_Cruncher_Konno.cfg'));

assert(triang.PARAM.SIG.FFTSIZE == 3840);
assert(isequal(triang.PARAM.SIG.T_LIM,[0 600]));
assert(isequal(triang.PARAM.SIG.F_LIM,[0.2 10]));
assert(triang.PARAM.SCRIPT.SMOOTHING_WIN_TYPE == 'T');
assert(konno.PARAM.SCRIPT.SMOOTHING_WIN_TYPE == 'K');
assert(konno.PARAM.KONNO.b == 40);

temporaryDirectory = tempname;
mkdir(temporaryDirectory);
cleanup = onCleanup(@() rmdir(temporaryDirectory,'s'));

validPath = fullfile(temporaryDirectory,'valid.cfg');
writeText(validPath,{ ...
    '% comment', ...
    'PARAM.A = -1.5e2; % inline comment', ...
    'PARAM.B = [1, 2; 3, 4];', ...
    'PARAM.C = ''value%retained'';', ...
    'PARAM.D = unquoted_text;'});
valid = readcfg(validPath);
assert(valid.PARAM.A == -150);
assert(isequal(valid.PARAM.B,[1 2; 3 4]));
assert(strcmp(valid.PARAM.C,'value%retained'));
assert(strcmp(valid.PARAM.D,'unquoted_text'));

assertInvalid({'PARAM.A = 1;','PARAM.A = 2;'}, ...
    temporaryDirectory,'duplicate.cfg');
assertInvalid({'PARAM.A = 1;','PARAM.A.B = 2;'}, ...
    temporaryDirectory,'value_then_group.cfg');
assertInvalid({'PARAM.A.B = 1;','PARAM.A = 2;'}, ...
    temporaryDirectory,'group_then_value.cfg');
assertInvalid({'PARAM.A [1,2];'}, ...
    temporaryDirectory,'missing_equals.cfg');
assertInvalid({'PARAM.A = [1,word];'}, ...
    temporaryDirectory,'invalid_numeric.cfg');
assertInvalid({'PARAM.A = ''unterminated;'}, ...
    temporaryDirectory,'unterminated_quote.cfg');
assertInvalid({'% comments only'}, ...
    temporaryDirectory,'empty.cfg','XNSR:ReadCfg:EmptyConfiguration');

fprintf('TUTTI I TEST READCFG SONO SUPERATI\n');

function assertInvalid(lines,directory,name,varargin)
path = fullfile(directory,name);
writeText(path,lines);
expectedIdentifier = 'XNSR:ReadCfg:InvalidLine';
if ~isempty(varargin)
    expectedIdentifier = varargin{1};
end
failedAsExpected = false;
try
    readcfg(path);
catch exception
    failedAsExpected = strcmp(exception.identifier,expectedIdentifier);
end
assert(failedAsExpected, ...
    'readcfg did not reject %s with the expected error.',name);
end

function writeText(path,lines)
fileId = fopen(path,'wt');
assert(fileId ~= -1,'Unable to create temporary CFG file.');
cleanup = onCleanup(@() fclose(fileId));
for lineIndex = 1:numel(lines)
    fprintf(fileId,'%s\n',lines{lineIndex});
end
clear cleanup
end
