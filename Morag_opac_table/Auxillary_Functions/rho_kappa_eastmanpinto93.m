function [ chi_1_exp , tau_nk] = rho_kappa_eastmanpinto93(nu , t, nu_nk , sig_nk,n_e,T)
% nu in units of energy
% n_e number density of electrons

sn_e = size(n_e);
if sn_e(2)>sn_e(1) % Prevents against 
    n_e=n_e';
end

if size(n_e,2)>1
    error('n_e is the wrong size')
end

c=set_consts();

tau_nk = c.c*c.h*t.*sig_nk.*n_e./nu_nk.*(1-exp(-nu_nk./T'));

[~, ~, groups] = histcounts(nu_nk, nu);
nu_grp_ind = groups > 1 & groups < (length(nu)-1);
zero_pad = zeros(length(nu)-max(groups)-1,1);

chi_1_exp = mid(nu)'./diff(nu)'/c.c/t.*[accumarray(groups(nu_grp_ind), 1-exp(-tau_nk(nu_grp_ind))) ; zero_pad];

chi_1_exp = chi_1_exp';
end