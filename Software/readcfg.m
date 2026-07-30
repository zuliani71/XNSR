function CFG = readcfg(cfgfile)
%READCFG Read an XNSR configuration file into a nested structure.
%   CFG = READCFG(FILE) accepts assignments such as:
%       PARAM.SIG.FFTSIZE = 3840;
%       PARAM.SIG.F_LIM = [0.2, 10];
%       PARAM.SCRIPT.SMOOTHING_WIN_TYPE = 'T';
%
%   Blank lines and MATLAB-style comments are ignored. Numeric scalars,
%   vectors and matrices are parsed without evaluating MATLAB expressions.

arguments
    cfgfile {mustBeTextScalar}
end

cfgfile = char(cfgfile);
assert(isfile(cfgfile), ...
    'XNSR:ReadCfg:MissingFile', ...
    'Configuration file not found: %s',cfgfile);

fileId = fopen(cfgfile,'rt');
assert(fileId ~= -1, ...
    'XNSR:ReadCfg:OpenFailed', ...
    'Unable to open configuration file: %s',cfgfile);
fileCleanup = onCleanup(@() fclose(fileId));
lines = textscan(fileId,'%s','Delimiter','\n','Whitespace','');
lines = lines{1};

CFG = struct();
seenFields = containers.Map('KeyType','char','ValueType','logical');

for lineNumber = 1:numel(lines)
    sourceLine = strtrim(lines{lineNumber});
    if isempty(sourceLine) || startsWith(sourceLine,'%')
        continue
    end

    sourceLine = strtrim(removeInlineComment(sourceLine));
    if isempty(sourceLine)
        continue
    end

    equalsIndex = find(sourceLine == '=',1,'first');
    if isempty(equalsIndex)
        cfgError(cfgfile,lineNumber,'Missing ''='' in assignment.');
    end

    fieldText = strtrim(sourceLine(1:equalsIndex-1));
    valueText = strtrim(sourceLine(equalsIndex+1:end));
    valueText = regexprep(valueText,';\s*$','');
    if isempty(fieldText) || isempty(valueText)
        cfgError(cfgfile,lineNumber,'Assignment requires a field and value.');
    end

    fields = strsplit(fieldText,'.');
    if any(cellfun(@(field) isempty(regexp(field, ...
            '^[A-Za-z][A-Za-z0-9_]*$','once')),fields))
        cfgError(cfgfile,lineNumber, ...
            sprintf('Invalid nested field name: %s',fieldText));
    end

    fieldKey = strjoin(fields,'.');
    if isKey(seenFields,fieldKey)
        cfgError(cfgfile,lineNumber, ...
            sprintf('Duplicate assignment for %s.',fieldKey));
    end

    try
        value = parseCfgValue(valueText);
    catch parseException
        if startsWith(parseException.identifier,'XNSR:ReadCfg:')
            cfgError(cfgfile,lineNumber,parseException.message);
        end
        rethrow(parseException);
    end

    CFG = add_nested_field(CFG,fields,value);
    seenFields(fieldKey) = true;
end

if seenFields.Count == 0
    error('XNSR:ReadCfg:EmptyConfiguration', ...
        'Configuration file contains no assignments: %s',cfgfile);
end

clear fileCleanup
end

function line = removeInlineComment(line)
inSingleQuote = false;
inDoubleQuote = false;
index = 1;
while index <= numel(line)
    character = line(index);
    if character == '''' && ~inDoubleQuote
        if inSingleQuote && index < numel(line) && line(index+1) == ''''
            index = index + 2;
            continue
        end
        inSingleQuote = ~inSingleQuote;
    elseif character == '"' && ~inSingleQuote
        if inDoubleQuote && index < numel(line) && line(index+1) == '"'
            index = index + 2;
            continue
        end
        inDoubleQuote = ~inDoubleQuote;
    elseif character == '%' && ~inSingleQuote && ~inDoubleQuote
        line = line(1:index-1);
        return
    end
    index = index + 1;
end
end

function value = parseCfgValue(valueText)
valueText = strtrim(valueText);
startsWithQuote = startsWith(valueText,'''') || startsWith(valueText,'"');
endsWithQuote = endsWith(valueText,'''') || endsWith(valueText,'"');
if xor(startsWithQuote,endsWithQuote) || ...
        (startsWithQuote && valueText(1) ~= valueText(end))
    error('XNSR:ReadCfg:InvalidValue', ...
        'Quoted text has unmatched delimiters.');
end

if numel(valueText) >= 2 && ...
        ((valueText(1) == '''' && valueText(end) == '''') || ...
         (valueText(1) == '"' && valueText(end) == '"'))
    quoteCharacter = valueText(1);
    value = valueText(2:end-1);
    value = strrep(value,[quoteCharacter quoteCharacter],quoteCharacter);
    return
end

if startsWith(valueText,'[') || endsWith(valueText,']')
    if ~(startsWith(valueText,'[') && endsWith(valueText,']'))
        error('XNSR:ReadCfg:InvalidValue', ...
            'Numeric array has unmatched brackets.');
    end
    value = parseNumericArray(valueText(2:end-1));
    return
end

if isNumericToken(valueText)
    value = numericTokenValue(valueText);
    return
end

if any(contains(valueText,{'[',']',',',';'}))
    error('XNSR:ReadCfg:InvalidValue', ...
        'Unsupported or malformed value: %s',valueText);
end

% Unquoted text is retained for compatibility with historical CFG files.
value = valueText;
end

function value = parseNumericArray(arrayText)
arrayText = strtrim(arrayText);
if isempty(arrayText)
    value = [];
    return
end

rowTexts = regexp(arrayText,';','split');
rows = cell(size(rowTexts));
columnCount = [];
for rowIndex = 1:numel(rowTexts)
    rowText = strtrim(rowTexts{rowIndex});
    if isempty(rowText)
        error('XNSR:ReadCfg:InvalidValue', ...
            'Numeric array contains an empty row.');
    end
    tokens = regexp(rowText,'[,\s]+','split');
    if any(cellfun(@(token) ~isNumericToken(token),tokens))
        error('XNSR:ReadCfg:InvalidValue', ...
            'Numeric array contains a nonnumeric token.');
    end
    rows{rowIndex} = cellfun(@numericTokenValue,tokens);
    if isempty(columnCount)
        columnCount = numel(rows{rowIndex});
    elseif numel(rows{rowIndex}) ~= columnCount
        error('XNSR:ReadCfg:InvalidValue', ...
            'Numeric array rows have inconsistent lengths.');
    end
end
value = vertcat(rows{:});
end

function valid = isNumericToken(token)
numericPattern = ['^[+-]?(?:(?:\d+\.?\d*)|(?:\.\d+))', ...
    '(?:[eEdD][+-]?\d+)?$'];
specialPattern = '^[+-]?(?:Inf|NaN)$';
valid = ~isempty(regexpi(token,numericPattern,'once')) || ...
    ~isempty(regexpi(token,specialPattern,'once'));
end

function value = numericTokenValue(token)
token = regexprep(token,'[dD]','e');
value = str2double(token);
end

function cfgError(cfgfile,lineNumber,message)
error('XNSR:ReadCfg:InvalidLine', ...
    '%s (file %s, line %d)',message,cfgfile,lineNumber);
end
