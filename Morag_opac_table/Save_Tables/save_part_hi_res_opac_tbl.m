%% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
%In Wexac it needs to be run in git
path_settings;

addpath(genpath(home_dir));
addpath(MG_opac_dir);
addpath([pathopac 'Kurucz/']);
addpath([pathopac 'nist_asd/']);

current_opac_table_settings;
% current_opac_table_test_settings;

% mixname = 'Solar'; Z_metal = [];
Xfrac_Sun = [0.7376 0.2494 0.0023865 0.0006902 0.005722 0.0012510 0.0007147 0.0006697 0.0003058 0.0012872];
Z_sun = sum(Xfrac_Sun(3:end));
mixname = 'Solar0_1Zs'; Z_metal = 0.1*Z_sun;
N_nu = 10^(6);
include_bf_bb = 1;
N_per_job = 4;

%% Consult the existing i_rho's
addpath(genpath(MG_opac_dir));
if exist([MG_opac_dir 'rhoT_exist_tbl.mat'],'file')
    load([MG_opac_dir 'rhoT_exist_tbl.mat'])
else
    rhoT_exist_tbl.started = zeros(1,length(tbl_rho)*length(tbl_T));
end
i_choice_start = find(~rhoT_exist_tbl.started,1,'first');
i_choice = i_choice_start + (0:N_per_job);
rhoT_exist_tbl.started(i_choice) = 1;
save([MG_opac_dir 'rhoT_exist_tbl.mat'],'rhoT_exist_tbl');

[i_rho, i_T] = ind2sub( [length(tbl_rho) length(tbl_T)] , i_choice);

%%
OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'i_start' num2str(i_choice_start) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];

%%
% Mode_temp.Plasma.A = A; % Atomic mass
% Mode_temp.Plasma.Z = Z; % Atomic number
% % Mode_temp.Plasma.Y = Y; % Degree of ionization (ne/np)
% Mode_temp.Plasma.Xfrac = Xfrac;
% Mode_temp.kappa.ff_on = 1;
% Mode_temp.kappa.bf_on = 1;
% Mode_temp.kappa.bb_on = 1;

if ~exist('OpacTableFilename','var')
    error('OpacTableFilename not defined')
end

tic
% profile on
[ kappa_abs,kappa_es,nu_calc ] = produce_part_high_res_opac_tbl( 1 , N_nu , tbl_rho , tbl_T ,include_bf_bb, Z_metal,i_rho , i_T); %Z_metal=[]->do nothing
% profile off
toc

% profsave(profile('info'),[profile_dir 'myprofile_produce_opac_table'])
save(OpacTableFilename , 'kappa_abs', 'kappa_es','nu_calc','-v7.3');

load([MG_opac_dir 'rhoT_exist_tbl.mat'])
rhoT_exist_tbl.finished(i_choice) = 1;
save([MG_opac_dir 'rhoT_exist_tbl.mat'],'rhoT_exist_tbl');