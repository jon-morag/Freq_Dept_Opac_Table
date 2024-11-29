
% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
%In Wexac it needs to be run in git
path_settings;

addpath(genpath(home_dir));
addpath(MG_opac_dir);
addpath([pathopac 'Kurucz/']);
addpath([pathopac 'nist_asd/']);

current_opac_table_RT_settings;
% current_opac_table_RT_settings_30_120;
% current_opac_table_test_settings;

% mixname = 'Solar'; Z_metal = [];
% Xfrac_Sun = [0.7376 0.2494 0.0023865 0.0006902 0.005722 0.0012510 0.0007147 0.0006697 0.0003058 0.0012872];
% Z_sun = sum(Xfrac_Sun(3:end));
% mixname = 'Solar0_1Z'; Z_metal = 0.1*Z_sun;
% mixname = 'Solar'; Z_metal = Z_sun;
% mixname = 'SolarZ0'; Z_metal = 0;
% 
% Xfrac_Sun = [0.7376 0.2494 0.0023865 0.0006902 0.005722 0.0012510 0.0007147 0.0006697 0.0003058 0.0012872];
% Z_sun = sum(Xfrac_Sun(3:end));
% % mixname = 'Solar0_1Z'; Z_metal = 0.1*Z_sun;
% mixname = 'Solar'; Z_metal = Z_sun;
% % mixname = 'SolarZ0'; Z_metal = 0;

% % C/O
% A = [12.0107 15.9994];
% nofrac = [0.5 0.5];
% Xfrac = nofrac.*A./(nofrac*A');
% Z = [6 8];
% mixname = 'CO';

% He-C/O
% A = [4.0026 12.0107 15.9994];
% COnofrac = [0.5 0.5];
% Z_CO = 0.1;
% Xfrac = [1-Z_CO Z_CO*COnofrac.*A(2:3)./(COnofrac*A(2:3)')];
% Z = [2 6 8];
% mixname = ['HeCO' num2str(Z_CO)];

mixname = 'H';

N_nu = 1e5;
low_res_spacing = 100;
include_bf = 0;
include_bb = 0;



% if include_bf_bb
%     OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];
% else
%     OpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff.mat'];
% end

OpacTableFilename = [MG_opac_dir mixname 'HiResOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) 'LR' num2str(low_res_spacing) '_es_ff'];
% OpacTableFilename = [MG_opac_dir 'no_HM_H_no_wings_' mixname 'HiResOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) 'LR' num2str(low_res_spacing) '_es_ff'];

if include_bf
    OpacTableFilename = [OpacTableFilename '_bf'];
end

if include_bb
    OpacTableFilename = [OpacTableFilename '_bb'];
end

Z_metal = 0;
% % OpacTableFilename = [OpacTableFilename '_rlines235' '.mat'];
% OpacTableFilename = [OpacTableFilename '_rlines5' '.mat'];

%
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

include_ff = 1;
tic
% profile on
[ kappa_abs,kappa_es,nu_calc ] = produce_high_res_opac_tbl_R_T_H_ff_only( 1 , N_nu , low_res_spacing, tbl_R , tbl_T,include_ff ,include_bf ,include_bb, Z_metal); %Z_metal=[]->do nothing
% [ kappa_abs,kappa_es,nu_calc ] = produce_high_res_opac_tbl_R_T_custom_mix( 1 , N_nu , low_res_spacing, tbl_R , tbl_T ,include_bf ,include_bb, A,Z,Xfrac); %Z_metal=[]->do nothing
% profile off
toc

% profsave(profile('info'),[profile_dir 'myprofile_produce_opac_table'])
if include_bf && include_bb
    save(OpacTableFilename , 'kappa_abs', 'kappa_es','nu_calc','tbl_T','tbl_R','-v7.3');
elseif include_bf
    kappa_abs_no_bb = kappa_abs;
    clear kappa_abs
    save(OpacTableFilename , 'kappa_abs_no_bb', 'kappa_es','nu_calc','tbl_T','tbl_R','-v7.3');
elseif include_bb
    error('error')
else
    kappa_abs_no_bf_bb = kappa_abs;
    clear kappa_abs
    save(OpacTableFilename , 'kappa_abs_no_bf_bb', 'kappa_es','nu_calc','tbl_T','tbl_R','-v7.3');
end