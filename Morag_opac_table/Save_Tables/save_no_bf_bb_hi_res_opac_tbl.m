%% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
%In Wexac it needs to be run in git
path_settings;

addpath(genpath(home_dir));
addpath(MG_opac_dir);
addpath([pathopac 'Kurucz/']);
addpath([pathopac 'nist_asd/']);

current_opac_table_settings;
% current_opac_table_test_settings;

mixname = 'Solar'; Z_metal = [];
% mixname = 'SolarZ0_1'; Z_metal = 0.1;
N_nu = 10^(5);
include_bf_bb = 0;

OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff.mat'];

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
[ kappa_abs_no_bf_bb,kappa_es,nu_calc ] = produce_high_res_opac_tbl( 1 , N_nu , tbl_rho , tbl_T ,include_bf_bb, Z_metal); %Z_metal=[]->do nothing

save(OpacTableFilename , 'kappa_abs', 'kappa_es','nu_calc');