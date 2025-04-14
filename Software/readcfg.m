function CFG=readcfg(cfgfile)
cfgData=textread(cfgfile,'%s','commentstyle','matlab','delimiter','\n');
unwantedChars = {';',''''};
CFG=[];
for i = 1:length(cfgData)
    lineCfg=char(cfgData{i});
    lineCfgSplit    = strsplit(lineCfg,'=');                % separate Parameters into name and value
    lineCfgVar      = strip(char(lineCfgSplit{1}));         % remove trailing blanks from Var name
    lineCfgVal      = strsplit(char(lineCfgSplit{2}),'%');  % separete Var value from possible comments
    lineCfgVal      = strip(char(lineCfgVal{1}));           % remove trailing blanks from Var value
    for j = 1:length(unwantedChars)
        lineCfgVal  = lineCfgVal(~ismember(lineCfgVal,...
            unwantedChars{j}));                             % remove unmwanted chars such;
    end
    lineCfgVarSplit = strsplit(lineCfgVar,'.');             % separate Var field names
    if ~isempty(str2num(lineCfgVal))
        Var=str2num(lineCfgVal);
    else
        Var= lineCfgVal;
    end
    CFG=add_nested_field(CFG,lineCfgVarSplit,Var);
end