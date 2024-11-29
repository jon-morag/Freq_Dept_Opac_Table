%% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
%In Wexac it needs to be run in git
path_settings;

addpath(genpath(home_dir));
addpath(MG_opac_dir);
addpath([pathopac 'Kurucz/']);
addpath([pathopac 'nist_asd/']);

current_opac_table_settings;
% current_opac_table_test_settings;

load('SolarMixOpacSpec2019.mat');
tbl_rho = SolarMixOpacSpec2019.rho;
tbl_T = SolarMixOpacSpec2019.T_eV;


% mixname = 'Solar'; Z_metal = [];
Xfrac_Sun = [0.7376 0.2494 0.0023865 0.0006902 0.005722 0.0012510 0.0007147 0.0006697 0.0003058 0.0012872];
Z_sun = sum(Xfrac_Sun(3:end));
% mixname = 'Solar0_1Zs'; Z_metal = 0.1*Z_sun;
mixname = 'Solar'; Z_metal = Z_sun;
include_bf_bb = 0;
low_res_spacing = 100;

% if include_bf_bb
%     OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];
% else
%     OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff.mat'];
% end

if include_bf_bb
    OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTableTOPSrho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(size(SolarMixOpacSpec2019.nu_eV,3)),2) '_es_ff_bf_bb.mat'];
else
    OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTableTOPSrho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(size(SolarMixOpacSpec2019.nu_eV,3)),2) '_es_ff.mat'];
end
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
% [ kappa_abs,kappa_es,nu_calc ] = produce_high_res_opac_tbl( 1 , N_nu , tbl_rho , tbl_T ,include_bf_bb, Z_metal); %Z_metal=[]->do nothing
% [ kappa_abs,kappa_es,nu_calc ] = produce_high_res_opac_tbl_R_rho( 1 , N_nu , tbl_rho , tbl_R ,include_bf_bb, Z_metal); %Z_metal=[]->do nothing
[ kappa_abs,kappa_es,nu_calc ] = produce_high_res_opac_tbl_TOPS_like( 1 ,low_res_spacing,SolarMixOpacSpec2019 ,include_bf_bb, Z_metal); %Z_metal=[]->do nothing
% profile off
toc


like_Solar_TOPS_2019 = SolarMixOpacSpec2019;
if include_bf_bb
    like_Solar_TOPS_2019.kappa_abs = kappa_abs;
else
    like_Solar_TOPS_2019.kappa_abs_no_bf_bb = kappa_abs;
    like_Solar_TOPS_2019 = rmfield(like_Solar_TOPS_2019,'kappa_abs');
end

% profsave(profile('info'),[profile_dir 'myprofile_produce_opac_table'])
save('like_Solar_TOPS_2019' , 'like_Solar_TOPS_2019' ,'-v7.3');
