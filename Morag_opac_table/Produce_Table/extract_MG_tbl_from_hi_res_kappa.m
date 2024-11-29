function [MG_tbl_out] = extract_MG_tbl_from_hi_res_kappa(Nnu_MG,nu_calc,kappa_abs,kappa_es,tbl_rho,tbl_T)

c=set_consts();
nu_lims_arr = logspace(-1,3,Nnu_MG)*c.eV; % erg
NG = Nnu_MG-1;

% Create multiple lower resolution groups from single opacity table.
kappa_Ross = zeros(length(tbl_T),length(tbl_rho),NG);
kappa_abs_avg = zeros(length(tbl_T),length(tbl_rho),NG);
kappa_planck_avg = zeros(length(tbl_T),length(tbl_rho),NG);

for i_rho = 1:length(tbl_rho)
    for i_T = 1:length(tbl_T)
        for i_nu_lim = 1:length(nu_lims_arr)-1
            seg_ind = nu_calc>=nu_lims_arr(i_nu_lim) & nu_calc<nu_lims_arr(i_nu_lim+1);
            nu_seg = nu_calc(seg_ind); nu_seg5 = mid(nu_seg);

            kappa_Ross(i_T,i_rho,i_nu_lim) = int_nu_dim(nu_seg,dBdT_EnergyDensity(tbl_T(i_T)*c.eV,nu_seg5),2) ./ int_nu_dim(nu_seg,dBdT_EnergyDensity(tbl_T(i_T)*c.eV,nu_seg5) ./ mid(reshape(kappa_abs(i_T,i_rho,seg_ind) + kappa_es(i_T,i_rho),[1,sum(seg_ind)])) , 2 );
            kappa_abs_avg(i_T,i_rho,i_nu_lim) = int_nu_dim(nu_seg, mid(reshape(kappa_abs(i_T,i_rho,seg_ind),[1,sum(seg_ind)])) , 2 )/(nu_lims_arr(i_nu_lim+1)-nu_lims_arr(i_nu_lim));
            kappa_planck_avg(i_T,i_rho,i_nu_lim) = int_nu_dim(nu_seg, PlanckEnergyDistribution(tbl_T(i_T)*c.eV,nu_seg5).*mid(reshape(kappa_abs(i_T,i_rho,seg_ind),[1,sum(seg_ind)])) , 2 ) ./ int_nu_dim(nu_seg, PlanckEnergyDistribution(tbl_T(i_T)*c.eV,nu_seg5) , 2 );
        end
    end
end

MG_tbl_out.kappa_Ross = kappa_Ross;
MG_tbl_out.kappa_avg = kappa_abs_avg;
MG_tbl_out.kappa_planck_avg = kappa_planck_avg;
MG_tbl_out.tbl_T = tbl_T;
MG_tbl_out.tbl_rho = tbl_rho;
MG_tbl_out.nu = mid(nu_lims_arr);