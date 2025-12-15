
function [ kappa_abs_hi_res , kappa_es, nu_calc ] = produce_hires_opac_tbl_rhoT(N_nu,tbl_rho ,tbl_T,include_ff,include_bf,include_bb,sobolev,expansion_line_limit,t_sobolev,A,Z,Xfrac )

c = set_consts();
nu_calc = logspace(-3.5,4,N_nu)*c.eV; % erg
% A later part of the calculation requires that nu_calc be log uniform

%% Extract Fields Mode Profile
Mode = standard_opac_profile();
Mode = load_required_tables(Mode);

if sobolev
    Mode.sobolev = 1;
    Mode.t_sobolev = t_sobolev; 
    Mode.nu_lims_arr = nu_calc; %average over groups so report the center
    nu_calc = mid(nu_calc);
    N_nu = N_nu-1;
end

if expansion_line_limit
    Mode.expansion_line_limit = 1;
    Mode.t_sobolev = t_sobolev;
end

Mode.Plasma.A = A; % Atomic mass
Mode.Plasma.Z = Z; % Atomic number
% Mode.Plasma.Y = Y; % Degree of ionization (ne/np)
Mode.Plasma.Xfrac = Xfrac;

if include_ff
    Mode.kappa.ff_on = 1;
else
    Mode.kappa.ff_on = 0;
end

if include_bf
    Mode.kappa.bf_on = 1;
else
    Mode.kappa.bf_on = 0;
end

if include_bb
    Mode.kappa.bb_on = 1;
else
    Mode.kappa.bb_on = 0;
end

if exist('Z_mix','var') && Z_mix == 0
    Mode.kappa.bb_Hydrogenic_species_select = 1:2;
end

if length(Mode.kappa.bf_species_select) >length(Z)
    Mode.kappa.bf_species_select = 1:length(Z);
    Mode.kappa.bb_species_select = 1:length(Z);
end

kappa_abs_hi_res = zeros(length(tbl_T),length(tbl_rho),length(nu_calc));
kappa_es = zeros(length(tbl_T),length(tbl_rho));

% table
for i_rho = 1:length(tbl_rho) % Run sequentially
% parfor i_rho = 1:length(tbl_rho) % Run in parallel
    for i_T = 1:length(tbl_T)
        if length(i_T)>1 && length(i_rho)>1
            disp(['Now computing rho ' num2str(i_rho) ' of ' num2str(length(tbl_rho)) ', T ' num2str(i_T) ' of ' num2str(length(tbl_T))])
        end
        [ne,Y,e_pop,Zeff] = SahaEq(tbl_T(i_T)*c.eV,tbl_rho(i_rho),Mode);
        if isnan(Y) && tbl_T(i_T) > 50 % Take care of high T limit
            Y=1-1e-3;
            for ix = 1:10
                e_pop(ix).z = [1e-7*ones(Mode.Plasma.Z(ix),1); 1-1e-3]; %non physical
            end
        end

        if isnan(Y) && tbl_T(i_T) < 1 % Take care of low T limit
            Y=1.9e-7;
            for ix = 1:10
                e_pop(ix).z = [1-1e-3; 1e-7*ones(Mode.Plasma.Z(ix),1)];
            end
        end
        kappa_abs_hi_res(i_T,i_rho,:)=kappa_abs_plasma(tbl_T(i_T)*c.eV,tbl_rho(i_rho),nu_calc,Mode,Y,e_pop);
        kappa_es(i_T,i_rho) = Y * sum(Mode.Plasma.Xfrac.*Mode.Plasma.Z./Mode.Plasma.A)*c.sigT/c.amu;
    end
end

if length(i_T)==1 && length(i_rho)==1
    kappa_abs_hi_res=squeeze(kappa_abs_hi_res);
    kappa_es = squeeze(kappa_es);
end

end