%% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
path_settings;
c = set_consts;

current_opac_table_RT_settings;
% current_opac_table_test_settings;
% mixname = 'Solar0_1Z';
% mixname = 'Solar';
% mixname = 'no_wings_Solar';
mixname = 'no_HM_H_no_wings_Solar';

% mixname = 'CO';
% Z_CO = 0.1; mixname = ['HeCO' num2str(Z_CO)];
% Z_CO = 0.5; mixname = ['HeCO' num2str(Z_CO)];


% mixname = 'SolarZ0';
N_nu = 1e6;
low_res_spacing = 100;
include_bf_bb = 1;

HiResOpacTableFilename = [MG_opac_dir mixname 'HiResOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) 'LR' num2str(low_res_spacing) '_es_ff_bf_bb.mat'];
load(HiResOpacTableFilename)

kappa_Ross = zeros(length(tbl_T),length(tbl_R));

Mode_temp.consts = set_consts();

for i_R = 1:length(tbl_R)
    for i_T = 1:length(tbl_T)
        kappa_Ross(i_T,i_R) = int_nu_dim(nu_calc,dBdT_EnergyDensity(tbl_T(i_T)*c.eV,mid(nu_calc),Mode_temp),2) ./ int_nu_dim(nu_calc,dBdT_EnergyDensity(tbl_T(i_T)*c.eV,mid(nu_calc),Mode_temp) ./ mid(reshape(kappa_abs(i_T,i_R,:) + kappa_es(i_T,i_R),[1,length(nu_calc)])) , 2 );
    end
end

OpacTableFilename = [MG_opac_dir  mixname 'RossMeanOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) 'LR' num2str(low_res_spacing) '_es_ff_bf_bb.mat'];
save(OpacTableFilename , 'kappa_Ross', 'kappa_es');