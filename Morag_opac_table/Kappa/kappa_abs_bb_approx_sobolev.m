function [k_bb] = kappa_abs_bb_approx_sobolev(nu,T,rho,Mode,F0,e_pop,Z,num_frac)
c=set_consts();
e_cutoff = 1e-4*c.eV; % Minimum transition energy

k_bb = zeros(1,length(nu));

for ix = 1:length(Z)
    if ismember(ix,Mode.kappa.bb_species_select)
        for ij = 1:Z(ix)
            ij_no_frac = e_pop(ix).z(ij)*num_frac(ix);
            if  ij_no_frac<Mode.Plasma.ion_pop_threshold || (Mode.Plasma.bb_i_state_select && Mode.Plasma.bb_i_state_select ~= ij)
                continue
            end

            if ((Z(ix)==1 && Mode.kappa.Hydrogen_bb_Nir) || (Z(ix)>1 &&  ij == Z(ix) && Mode.kappa.Hydrogenic_bb_Nir)) && ismember(ix,Mode.kappa.bb_Hydrogenic_species_select) %Hydrogen or Hydrogenic!
                [n_n,nunk,DecayRate,signk]  = signk_bb_ij_Hlike(T,Mode,Z(ix),F0,c,e_cutoff);
            elseif ij ~= Z(ix) && Mode.kappa.bb_Kurucz && ~isempty(Mode.Transition_Rate_Matrix(Z(ix),ij).fnk)
                [n_n,nunk,DecayRate,signk] = signk_bb_ij_Kurucz(ij,T,Mode,Z(ix),F0,c,e_cutoff);
            else
                continue
            end
            
%             % Find line index, assuming nu is logarithmically uniform
%             nunk_ind = round((length(nu_hi_res)-1)/(l_max_nu-l_min_nu)*(log(nunk)-l_min_nu) + 1);

%             for i_line=1:NL
%                 if ij_no_frac*n_n(i_line)<Mode.Plasma.ion_pop_threshold
%                     continue
%                 end
% 
%                 if nunk_ind(i_line) < 1 || nunk_ind(i_line) > length(nu_hi_res)
%                     continue
%                 end
% 
%                 if nunk(i_line) < nu(1) || nunk(i_line) > nu(end)
%                     continue
%                 end
% 
%             end

                [chi_1 tau_s] = rho_kappa_eastmanpinto93(Mode.nu_lims_arr , Mode.t_sobolev, nunk , signk, rho.*ij_no_frac*n_n/c.amu,T);
                k_bb = k_bb + chi_1/rho;

            disp(['Computed bb for Z ', num2str(Z(ix)), ' ij ', num2str(ij)] )

        end
    end
end
end

function [n_n, nu_nk,DecayRate,signk] = signk_bb_ij_Hlike(T,Mode,Z,F0,c,e_cutoff)
IH = c.IH;
nlevels_Hbb=min(length(Mode.H_Transition_Rate_Matrix.fnk),Mode.Ionization.nmax_levels_Hydrogenic);
fnk = Mode.H_Transition_Rate_Matrix.fnk(1:nlevels_Hbb,1:nlevels_Hbb);
nu_nk = Z^2*Mode.H_Transition_Rate_Matrix.nu_nk(1:nlevels_Hbb,1:nlevels_Hbb);
DecayRate = Z^4*Mode.H_Transition_Rate_Matrix.DecayRate(1:nlevels_Hbb,1:nlevels_Hbb); %% I removed the A(1)/A(ix) dependence.
ExRateCoef = 0; %Currently not implemented.

%                     %% For looking at the population levels
%                     figure(1)
%                     semilogy(squeeze(n_nM(1,:,1)),'DisplayName',['T = ' num2str(Tm5n0/c.eV) 'eV, \rho = ' num2str(rhom5n0) 'gcm^{-3}'])
%                     hold on
%                     xlabel('quantum number n'); ylabel('Population fraction');
%                     ax = gca(); ax.FontSize=14;
%                     legend show

indL = nu_nk >= e_cutoff; 

nhi = ones(nlevels_Hbb,1) * (1:nlevels_Hbb) ;
nlo = nhi';
nhi = nhi(indL);
nlo = nlo(indL);
fnk = fnk(indL);
nu_nk = nu_nk(indL);
DecayRate = DecayRate(indL);

gn_lo = 2*nlo.^2;
En = IH*Z^2./nlo.^2;
n_sq = (1:nlevels_Hbb).^2;

if Mode.kappa.include_w_hm_Hydrogenic
    w_hm = get_w_HM(F0,nlo,1,1,Mode,c,Z);
    Zr = gn_lo.*w_hm.*exp(-(IH*Z^2-En)./T');
    w_hm = get_w_HM(F0,1:nlevels_Hbb,1,1,Mode,c,Z); %Line Index is different so must repeat w_hm.
    Gr = sum(2*n_sq.*w_hm.*exp(-IH*Z^2*(1-1./n_sq)./T'),2);
else
    Zr = gn_lo.*exp(-(IH*Z^2-En)./T');
    Gr = sum(2*n_sq.*exp(-IH*Z^2*(1-1./n_sq)./T'),2);
end

n_n = Zr./Gr;
signk = pi*c.q^2/c.me/c.c*fnk;

end

function [n_n, nu_nk, DecayRate,signk] = signk_bb_ij_Kurucz(ij,T,Mode,Z,F0,c,e_cutoff)
%% Get line info from matrix
nlines=min(ceil(length(Mode.Transition_Rate_Matrix(Z,ij).fnk)*Mode.kappa.levels_cutoff_frac),Mode.kappa.nmax_lines_Kurucz);
% levels cutoff frac is to test our sensitivity to the finite number of
% levels

fnk = Mode.Transition_Rate_Matrix(Z,ij).fnk(1:nlines);
nu_nk = Mode.Transition_Rate_Matrix(Z,ij).nu_nk(1:nlines);
DecayRate = Mode.Transition_Rate_Matrix(Z,ij).DecayRate(1:nlines);
%                         ExRateCoef = Mode.Transition_Rate_Matrix(Z(ix),ij).ExRateCoef(1:nlines)';
ExRateCoef = 0; %Currently not supported
gn_lo = Mode.Transition_Rate_Matrix(Z,ij).glo(1:nlines);
if isfield(Mode.Transition_Rate_Matrix(Z,ij),'ghi')
    gn_hi = Mode.Transition_Rate_Matrix(Z,ij).ghi(1:nlines);
end
E_n = Mode.Transition_Rate_Matrix(Z,ij).Elo(1:nlines)';
%Should not use Enhi. Way too many NaN's in Kurucz. Extract
%from Enlo instead.

ind_e = nu_nk >= e_cutoff; % need to check this cutoff vs. opacity of highly
%excited lines, which could be low due to colliosional broadening

fnk = fnk(ind_e);
nu_nk = nu_nk(ind_e);
DecayRate = DecayRate(ind_e);
gn_lo = gn_lo(ind_e);
E_n = E_n(ind_e);
I = Mode.Ionization.energies(Z,ij)*c.eV;

%% Build Partition function
if Mode.kappa.include_w_hm_nonHydrogenic
    n_eff = ij*real(sqrt(IH./(I-E_n)));
    w_hm = get_w_HM(F0,n_eff,1,1,Mode,c,ij);
    w_hm(n_eff==0) = 0;
else
    w_hm = 1;
end

Zr = gn_lo'.*w_hm.*exp(-E_n./T');

% For low temperatures:
% if any(any(Zr==0)); Zr = gn_lo'.*exp(-(I+E_n)./Tm5n0' + min((I+E_n)./Tm5n0',[],2)); end

% Partition function only over unique beginning states
v = [gn_lo, E_n'/c.eV];
v = sort(v, 2);
[uniquerows, idx] = unique(v, 'rows');
Gr = sum(Zr(:,idx),2);

%% Optional Artificially extend the partition function to check our sensitivity to finite energy levels
if Mode.kappa.extend_partition_func_n
    n_extend = round(max(n_eff)) + 1 : Mode.kappa.extension_quantum_number;
    En_extend = I - ij^2*c.IH./n_extend.^2;

    if Mode.kappa.include_w_hm_nonHydrogenic
        w_hm_extend = get_w_HM(F0,n_extend,1,1,Mode,c,ij);
        Gr_extend = sum(2*n_extend.^2.*w_hm_extend.*exp(-En_extend./T'));
    else
        Gr_extend = sum(2*n_extend.^2.*exp(-En_extend./T'));
    end
    Gr = Gr + Gr_extend;
end

%% Finalize
n_n = Zr./Gr;
signk = pi*c.q^2/c.me/c.c*fnk;

if Mode.kappa.fix_A_decay_Kurucz && exist('gn_hi')
    %Kurucz doesn't always provide Decay rates, so assume Ank as minimum

    Ank = 8*pi*gn_lo./gn_hi.*(nu_nk/c.h).^2.*signk/c.c^2;
    DecayRate(DecayRate==0) = Ank(DecayRate==0);
end

end

function [out] = diff2(X)
    out = [X(2)-X(1),(X(3:end)-X(1:end-2))/2,X(end)-X(end-1)];
end