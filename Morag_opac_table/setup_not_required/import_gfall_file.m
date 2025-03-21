%% Import data from text file.
% Script for importing data from the following text file:
%
%    D:\Dropbox (Weizmann Institute)\OpacityTables\gfall.dat
%
% To extend the code to different selected data or a different text file, generate a function instead of a script.

% Auto-generated by MATLAB on 2023/06/09 08:40:33

function gfall_tbl = import_gfall_file(save_dir,filename_import ,filename_save)
%% Initialize variables.
startRow = 3;

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: categorical (%C)
%   column7: double (%f)
%	column8: double (%f)
%   column9: categorical (%C)
%	column10: double (%f)
%   column11: double (%f)
% For more information, see the TEXTSCAN documentation.
% formatSpec = '%11f%7f%6f%12f%6f%11C%11f%5f%11C%6f%6f%[^\n\r]';
  formatSpec = '%11f%7f%6f%12f%6f%10C%12f%5f%11C%6f%6f%[^\n\r]';

% only import the 1st column

%% Open the text file.
fileID = fopen([save_dir filename_import],'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

%% Create output variable
gfall_tbl = table(dataArray{1:end-1}, 'VariableNames', {'wl_nm','log_GF','elem_chrg','E1_cm_1','J1','label1','E2_cm_1','J2','label2','log_Gam_rad','log_Gam_Strk'});
save([save_dir filename_save],'gfall_tbl','-v7.3');
end