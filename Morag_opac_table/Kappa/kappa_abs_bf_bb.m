function [k_bf_bb]=kappa_abs_bf_bb(T,rho,nu,Mode,Y,e_pop)
% calculates bf and bb for atomic mix
% bf: Analytic for H-like (one electron), Verner Yakovlev table otherwise
% bb: Analytic for H-like, Kurucz otherwise
% for efficiency, bb line peaks are computed in hi-res, and line wings are
% computed with lower freq. res 

% Assumes that nu is logarithmically uniform. Will not work correctly otherwise!!

%% Setup
c = Mode.consts;

% nmax = min(max(nmaxm5n0),10);
if ~isfield(Mode.kappa,'nmax_lines_Kurucz')
    Mode.kappa.nmax_lines_Kurucz = 1000; % default value
end

% Most calculations will be done in low resolution, except the center of
% lines:
nu_hi_res = nu;
low_res_spacing = Mode.kappa.low_res_spacing; % How low the resolution is relative to 
low_res_ind = [1:low_res_spacing:(length(nu_hi_res)-1) length(nu_hi_res)];
nu_low_res = nu(low_res_ind);
nu = nu_low_res;
 
n = 1:Mode.Ionization.nmax_levels_Hydrogenic;

X = Mode.Plasma.Xfrac;
Z = Mode.Plasma.Z;
A = Mode.Plasma.A;

natomm5n0 = rho*sum(X./A)/c.amu;
rion = (4*pi*natomm5n0/3).^(-1/3);
F0 = Y*sum(X.*Z./A)/sum(X./A)./rion.^2;

num_frac = X./A/sum(X./A);

%% Bound-free
if Mode.kappa.bf_on
    k_bf = kappa_abs_bf(n, T, nu, Mode,e_pop,Z,num_frac,F0,c);
else
    k_bf = 0;
end

%% Bound-bound

if Mode.kappa.bb_on
    [k_bb, k_bb_hi_res] = kappa_abs_bb(nu_hi_res,nu,T,Mode,F0,e_pop, Z,A, X./A/sum(X./A));
else
    k_bb = 0;
    k_bb_hi_res = 0;
end

%% Combine
k_bb_hi_res(k_bb_hi_res<0) = 0; %Due to the subtraction we perform, there can be a negative value due to poor interp.
k_bf_bb = interp1( nu , k_bf + k_bb , nu_hi_res) + k_bb_hi_res;

end