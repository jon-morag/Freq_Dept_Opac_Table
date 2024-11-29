%% Setup
req_tbls_dir = fullfile([fileparts(mfilename('fullpath')) '/../Required_Tables/']);
gfall_filename = 'gfall';
max_element = 26; % i.e. we choose to stop extracting the table at Fe
matrix_save_name = 'Kurucz_lines_cd23';

%% Choose Kurucz line list
% CD 1 not supported
% CD 23
line_url = 'http://kurucz.harvard.edu/linelists/gfall/gfall.dat';
% % Most recent url: (also not supported, user should import manually)
% line_url = 'http://kurucz.harvard.edu/linelists/gfnew/gffall08oct17.dat';

%% Produce Transition Rate Matrix
disp('Downloading from Kurucz')
websave([req_tbls_dir gfall_filename '.dat'],line_url);
disp('Putting Kurucz list in Matrix Format')
bb_lines_tbl = import_gfall_file(req_tbls_dir,[gfall_filename '.dat'] ,[gfall_filename '_tbl.dat']);
Transition_Rate_Matrix = produce_save_bb_mat_from_Kurucz(bb_lines_tbl,max_element,[req_tbls_dir matrix_save_name '.mat']);

