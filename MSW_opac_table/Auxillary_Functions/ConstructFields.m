var_names=who;
% See also ExtractFields
clear All_vars;
for l=1:length(var_names);
eval(['All_vars.',var_names{l},'=',var_names{l},';']);
end;
