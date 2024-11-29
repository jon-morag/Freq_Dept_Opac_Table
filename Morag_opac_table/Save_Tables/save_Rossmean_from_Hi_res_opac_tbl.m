%% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
path_settings;
c = set_consts;

current_opac_table_settings;
% current_opac_table_test_settings;
% mixname = 'Solar0_1Z';
mixname = 'Solar';
N_nu = 10^(5);

HiResOpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];
load(HiResOpacTableFilename)

kappa_Ross = zeros(length(tbl_T),length(tbl_rho));

Mode_temp.consts = set_consts();

for i_rho = 1:length(tbl_rho)
    for i_T = 1:length(tbl_T)
        kappa_Ross(i_T,i_rho) = int_nu_dim(nu_calc,dBdT_EnergyDensity(tbl_T(i_T)*c.eV,mid(nu_calc),Mode_temp),2) ./ int_nu_dim(nu_calc,dBdT_EnergyDensity(tbl_T(i_T)*c.eV,mid(nu_calc),Mode_temp) ./ mid(reshape(kappa_abs(i_T,i_rho,:) + kappa_es(i_T,i_rho),[1,length(nu_calc)])) , 2 );
    end
end

OpacTableFilename = [MG_opac_dir  mixname 'RossMeanOpacTablerho' num2str(length(tbl_rho)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];
save(OpacTableFilename , 'kappa_Ross', 'kappa_es');