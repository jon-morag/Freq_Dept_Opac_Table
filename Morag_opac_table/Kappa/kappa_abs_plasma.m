function [k_abs]=kappa_abs_plasma(T,rho,nu,Mode,Ym5n0,e_pop)
% opacity for an arbitrary plasma with Y degree of ionization - including
% free-free, free-bound and bound-free emission. nu is in ergs


if Mode.kappa.ff_on
    k_ff = kappa_abs_Brems(T,rho,nu,Mode,Ym5n0,e_pop);
else
    k_ff = 0;
end

if isfield(Mode,'sobolev') && Mode.sobolev
    k_bf_bb = kappa_abs_bf_bb_approx_sobolev(T,rho,nu,Mode,Ym5n0,e_pop);
else
    k_bf_bb = kappa_abs_bf_bb(T,rho,nu,Mode,Ym5n0,e_pop);
    % k_bf_bb = kappa_abs_bf_bb_old_prior_phi_func(T,rho,nu,Mode,Ym5n0,e_pop); %test
end

k_abs = k_ff+k_bf_bb;
