function [k_bf_bb]=kappa_abs_bf_bb_approx_sobolev(T,rho,nu,Mode,Y,e_pop)
% calculates bf and bb for atomic mix
% bf: Analytic for H-like (one electron), Verner Yakovlev table otherwise
% bb: Sobolev (Eastman Pinto 93)
% for efficiency, bb line peaks are computed in hi-res, and line wings are
% computed with lower freq. res 

% Assumes that nu is logarithmically uniform. Will not work correctly otherwise!!

%% Setup
c = Mode.consts;

% nmax = min(max(nmaxm5n0),10);
if ~isfield(Mode.kappa,'nmax_lines_Kurucz')
    Mode.kappa.nmax_lines_Kurucz = 1000; % default value
end

% Unlike the high-res case, all calculations will be done in low res nu, so
% no need to split nu into hi- and low-res.

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
    [k_bb] = kappa_abs_bb_approx_sobolev(nu,T,rho,Mode,F0,e_pop, Z, X./A/sum(X./A));
else
    k_bb = 0;
end

%% Combine
k_bf_bb = k_bf + k_bb;

end