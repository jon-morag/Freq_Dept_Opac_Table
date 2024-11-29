function [u_low_res_mapped,nu_low_res_edge,nu_hi_res_sample,u_hi_res_mapped,len_phi,indcs,nu_mid_hi_res_sample,u_mid_hi_res,u_mid_hi_res_mapped]= split_low_hi_res_map_to_u_fconsv(nu,nu_hi_res,nu_mid_hi_res,nunk,Dnu,n_x_FWHM,n_x_FWHM_fconsv,min_u,max_u,len_u,low_res_spacing)
u_low_res_mapped=[];
nu_low_res_edge=[];
nu_hi_res_sample=[];
u_hi_res_mapped=[];
len_phi=[];
indcs=[];

nu_mid_hi_res_sample=[];
u_mid_hi_res=[];
u_mid_hi_res_mapped=[];

nu_low_edge_idx = find(nu > nunk-max_u/Dnu,1,'first');
nu_low_bnd_edge_idx = find(nu(nu_low_edge_idx:end) < (nunk -  n_x_FWHM),1,'last');
nu_low_bnd_idx = nu_low_edge_idx + nu_low_bnd_edge_idx-1;
nu_hi_bnd_idx = nu_low_bnd_idx + find(nu(nu_low_bnd_idx+1:end) > (nunk +  n_x_FWHM),1,'first');
if isempty(nu_low_bnd_idx) || isempty(nu_hi_bnd_idx)
    return
end
nu_hi_edge_idx = nu_hi_bnd_idx + find(nu(nu_hi_bnd_idx+1:end) < nunk+max_u/Dnu,1,'last');
nu_hi_bnd_edge_idx = nu_hi_bnd_idx - nu_low_edge_idx+1;
len_phi = 1+nu_hi_edge_idx-nu_low_edge_idx;

low_res_bin_idx = [1:(1+nu_low_bnd_idx-nu_low_edge_idx) (1+nu_hi_bnd_idx-nu_low_edge_idx):len_phi];
nu_low_res_edge = nu(nu_low_edge_idx:nu_hi_edge_idx);
nu_low_res_sample = nu([nu_low_edge_idx:nu_low_bnd_idx nu_hi_bnd_idx:nu_hi_edge_idx]);

nu_hi_res_low_idx = 2 + (nu_low_bnd_idx - 1)*low_res_spacing;
nu_hi_res_hi_idx = (nu_hi_bnd_idx - 1)*low_res_spacing;
hi_res_bin_idx = nu_hi_res_low_idx : nu_hi_res_hi_idx;

nu_hi_res_sample = nu_hi_res(hi_res_bin_idx);
% fconsv only
f_cons_idx = abs(nu_hi_res_sample-nunk) <  n_x_FWHM_fconsv;
nu_hi_res_sample_f_cons = nu_hi_res_sample(~f_cons_idx);
%

u_low_res = (nu_low_res_sample-nunk).*Dnu;
u_low_res(abs(u_low_res)<min_u) = min_u; % We don't care on which side of zero it is. Not sensitive.
u_low_res_mapped = (len_u-1)/(log(max_u)-log(min_u))*(log(abs(u_low_res))-log(min_u)) + 1;

u_hi_res = (nu_hi_res_sample_f_cons-nunk).*Dnu;
u_hi_res_mapped = (len_u-1)/(log(max_u)-log(min_u))*(log(abs(u_hi_res))-log(min_u)) + 1;
% u_hi_res_mapped(u_hi_res_mapped>len_u) = len_u; not in fconst
u_hi_res_mapped(u_hi_res_mapped<1) = 1;

% fconsv only
hi_res_mid_bin_idx = (nu_hi_res_low_idx-1) : nu_hi_res_hi_idx;

nu_mid_hi_res_sample = nu_mid_hi_res(hi_res_mid_bin_idx);
nu_mid_hi_res_sample = nu_mid_hi_res_sample(find(f_cons_idx + [0 f_cons_idx(1:end-1)]));
if f_cons_idx(end)
    nu_mid_hi_res_sample = [nu_mid_hi_res_sample nu_mid_hi_res(nu_hi_res_hi_idx)];
end

u_mid_hi_res = (nu_mid_hi_res_sample-nunk).*Dnu;
u_mid_hi_res_mapped = (len_u-1)/(log(max_u)-log(min_u))*(log(abs(u_mid_hi_res))-log(min_u)) + 1;
u_mid_hi_res_mapped(u_mid_hi_res_mapped<1) = 1;

indcs.low_res_bin = low_res_bin_idx;
indcs.hi_res_bin = hi_res_bin_idx;
indcs.nu_low_bnd = nu_low_bnd_idx;
indcs.nu_low_bnd_edge = nu_low_bnd_edge_idx;
indcs.nu_low_edge = nu_low_edge_idx;
indcs.nu_hi_edge = nu_hi_edge_idx;
indcs.nu_hi_bnd = nu_hi_bnd_idx;
indcs.nu_hi_bnd_edge = nu_hi_bnd_edge_idx;

indcs.f_cons=f_cons_idx;

end