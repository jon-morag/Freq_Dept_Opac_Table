%% Is now free standing. Can be used without getting Mode from PolytropeDiffusionMGSahaMix
path_settings;
c = set_consts;

current_opac_table_RT_settings;
% current_opac_table_test_settings;
% mixname = 'Solar0_1Z';
mixname = 'Solar';
N_nu = 10^(5);

HiResOpacTableFilename = [MG_opac_dir  mixname 'HiResOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];
load(HiResOpacTableFilename)

kappa_avg = zeros(length(tbl_T),length(tbl_R));

Mode_temp.consts = set_consts();

for i_R = 1:length(tbl_R)
    for i_T = 1:length(tbl_T)
        kappa_avg(i_T,i_R) = int_nu_dim(nu_calc,PlanckEnergyDensity(tbl_T(i_T)*c.eV,mid(nu_calc),Mode_temp).* mid(reshape(kappa_abs(i_T,i_R,:),[1,length(nu_calc)])), 2 )./int_nu_dim(nu_calc,PlanckEnergyDensity(tbl_T(i_T)*c.eV,mid(nu_calc),Mode_temp),2);
    end
end

OpacTableFilename = [MG_opac_dir  mixname 'BBAvgOpacTableR' num2str(length(tbl_R)) 'T' num2str(length(tbl_T)) 'Nnu1e' num2str(log10(N_nu),2) '_es_ff_bf_bb.mat'];
save(OpacTableFilename , 'kappa_avg');