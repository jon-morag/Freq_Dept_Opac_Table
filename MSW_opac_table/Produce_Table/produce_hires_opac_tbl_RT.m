
function [ kappa_abs_hi_res , kappa_es, nu_calc ] = produce_hires_opac_tbl_RT(N_nu,tbl_R ,tbl_T,include_ff,include_bf,include_bb,A,Z,Xfrac )

c = set_consts();
nu_calc = logspace(-3.5,4,N_nu)*c.eV; % erg

kappa_abs_hi_res = zeros(length(tbl_T),length(tbl_R),length(nu_calc));
kappa_es = zeros(length(tbl_T),length(tbl_R));

%% Extract Fields Mode Profile
Mode = standard_opac_profile();
Mode = load_required_tables(Mode);

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

% table
for i_R = 1:length(tbl_R) % Run sequentially
% parfor i_R = 1:length(tbl_R) % Run in parallel
    for i_T = 1:length(tbl_T)
        disp(['Now computing R ' num2str(i_R) ' of ' num2str(length(tbl_R)) ', T ' num2str(i_T) ' of ' num2str(length(tbl_T))])
        tbl_rho = ( tbl_T(i_T) * c.eV ).^3./tbl_R(i_R);
        [ne,Y,e_pop,Zeff] = Saha(tbl_T(i_T)*c.eV,tbl_rho,Mode);
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
        kappa_abs_hi_res(i_T,i_R,:)=kappa_abs_plasma(tbl_T(i_T)*c.eV,tbl_rho,nu_calc,Mode,Y,e_pop);
        kappa_es(i_T,i_R) = Y * sum(Mode.Plasma.Xfrac.*Mode.Plasma.Z./Mode.Plasma.A)*c.sigT/c.amu;
    end
end

end