function [k_shifted] = doppler_broaden_shift(nu,kappa,v_doppler,Deltav_doppler)
% Approximate a doppler broadening and shifting of a hi-res opacity array or
% an entire hi-res table
% assumes nu is loguniform!

% careful of the RAM requirement in large opacity tables (e.g. Nnu>1e6)

%Input v_doppler (shift), Deltav_doppler (broaden) in cm/s
c=set_consts();

dnu_nu = v_doppler/c.c;
Deltanu_nu = Deltav_doppler/c.c;
alpha_nu = nu(2)/nu(1);

k_shifted = fracbroadenshift(kappa,log(1-Deltanu_nu)/log(alpha_nu) , log(1+dnu_nu)/log(alpha_nu));

end