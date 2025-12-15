function [k_bb,k_bb_hi_res] = kappa_abs_bb(nu_hi_res,nu,T,Mode,F0,e_pop,Z,A,num_frac)
c=set_consts();
e_cutoff = 1e-4*c.eV; % Minimum transition energy

diff_nu = diff2(nu_hi_res);
nu_mid_hi_res = [nu_hi_res(1) , (nu_hi_res(2:end) + nu_hi_res(:,1:end-1,:))/2 , nu_hi_res(end)];
H_table = Mode.HjertingTable;
l_max_a = log(max(H_table.a));
l_min_a = log(min(H_table.a));
max_u = max(H_table.u);
min_u = min(H_table.u);
l_max_nu = log(max(nu));
l_min_nu = log(min(nu));

k_bb = zeros(1,length(nu));
k_bb_hi_res = zeros(1,length(nu_hi_res));

if isfield(Mode,'expansion_line_limit') && Mode.expansion_line_limit
    % limit only the hi-res part (peak not line wings) to be expansion
    % limited. Limit the n sigma part but do the nu multiplication later
    kappa_l_sob_lim = Mode.kappa_l_sob_lim;
    expansion_limit = 1;
else
    expansion_limit = 0;
end


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

            NL = length(nunk); % number of lines


            %% Physics Extensions
            % Electron Impact Broadening

            if Z(ix) == ij && Mode.kappa.H_e_impact_broadening %H-like only!
                %Full half - width %Armstrong 1964
                nem5n0 = Y.*rho*(Z./A*X')/c.amu;
                FWHM_cl = 8*pi^2/Z(ix)^2*c.q^2*c.ale*c.a0^2*nem5n0.*sqrt(c.me*c.c^2./(2*pi*T)).* nhiL.^4;
                %                     n_critical = (sqrt(Tm5n0/32/pi^3/c.me/c.c/c.c)/(c.ale*c.a0^3*nem5n0)).^(1/7);
                %                     DnuDMm5n0_1 = ( 1./DnuDMm5n0_1.^2 + FWHM_cl.^2/4/log(2) ).^(-1/2);

                DecayRate = DecayRate + 2*pi*FWHM_cl/c.h;
                %                     DecayRate = sqrt(DecayRate.^2 + (2*pi*FWHM_cl/c.h).^2);
            end

            if ij == Z(ix) && Mode.kappa.Hydrogenic_Dappen_HM_factor
                w_hm_DAM_bb = get_w_HM(F0,nhiL,1,1,Mode,c,ij);
            elseif Z(ix) > ij && Mode.kappa.Kurucz_Dappen_HM_factor
                nhi_eff = ij*real(sqrt(IH./(I-E_n(indL)-nunk')));
                % Some of the E_nhi's are not a number. Use E_n and
                % nunk instead
                w_hm_DAM_bb = get_w_HM(F0,nhi_eff,1,1,Mode,c,ij);
                w_hm_DAM_bb(nhi_eff==0) = 0;
            else
                w_hm_DAM_bb = ones(1,NL);
            end

            %% FWHM for the Voigt Function
            vth = (2*T/A(ix)/c.amu).^(1/2);
            Dnu = c.c./(nunk*vth);

            ExRate=0; % Currently not supported
            FWHM_L = (DecayRate+2*ExRate)/2/pi;
            FWHM_G = sqrt(log(2))/c.h./Dnu;
            FWHM_voigt = FWHM_L/2 + sqrt(FWHM_L.^2/4 + FWHM_G.^2); %approximate formula for Voigt
            
            % Holtsmark distribution 'a' parameter:
            a = c.h/2*FWHM_L.*Dnu;
            a(a<min(H_table.a)) = min(H_table.a);
            %                 a(a>max(H_table.a)) = max(H_table.a);
            a_mapped = (length(H_table.a)-1)/(l_max_a-l_min_a)*(log(a)-l_min_a) + 1;

            if isfield(Mode.kappa,'plot_FWHM') && Mode.kappa.plot_FWHM
                p_Z(plot_i) = semilogy(FWHM_L,'DisplayName',['Z = ' num2str(Z(ix)) ', Ionization =' num2str(ij)]);
                hold on
                ax=gca(); ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
                semilogy(FWHM_G,'LineWidth',3);
                [dummy nunk_ind] = min(abs(nunk-nu)');
                ax=gca(); ax.ColorOrderIndex=mod(ax.ColorOrderIndex-2,7)+1;
                semilogy(nu(nunk_ind)/c.h,'--','LineWidth',2);

                p_L = plot([NaN, NaN],'k','DisplayName','Lorenzian'); p_G = plot([NaN, NaN],'k','LineWidth',3,'DisplayName','Gaussian'); p_dnu = plot([NaN, NaN],'k--','LineWidth',2,'DisplayName','\Delta \nu (grid)');
                legend([p_Z p_L p_G p_dnu]);

                title(['\rho=', num2str(rho(1)) 'g cm^{-3}, T=' num2str(T/c.eV) 'eV'])
                xlabel('line Number');ylabel('FWHM [s^{-1}]');
                plot_i = plot_i+1;
            end

            
           %% Produce voight function for each line
           low_res_spacing = Mode.kappa.low_res_spacing; % How low the resolution is relative to 
           n_FWHM_hi_res = Mode.kappa.n_FWHM_hi_res; % Number of FWHM's from the peak to calc in hi-res
            % Find line index, assuming nu is logarithmically uniform
            nunk_ind = round((length(nu_hi_res)-1)/(l_max_nu-l_min_nu)*(log(nunk)-l_min_nu) + 1);

            for i_line=1:NL
                if ij_no_frac*n_n(i_line)<Mode.Plasma.ion_pop_threshold
                    continue
                end

                if nunk_ind(i_line) < 1 || nunk_ind(i_line) > length(nu_hi_res)
                    continue
                end

                if strcmp(Mode.kappa.binning_function,'Optimal')
                    % Shape determined by line width relative to grid
                    if FWHM_voigt(i_line)<diff_nu(nunk_ind(i_line))/c.h/10 && FWHM_L(i_line)<diff_nu(nunk_ind(i_line))/c.h/3e3 && Z(ix) > ij% Keeps certain lines with wide wings that may affect the continuum.
                        binning_function = 'DeltaOnly';
                    elseif FWHM_voigt(i_line) > diff_nu(nunk_ind(i_line))/c.h*3
                        binning_function = 'SampleAtBinEdge';
                    else
                        binning_function = 'fConserving';
                    end
                else
                    binning_function = Mode.kappa.binning_function;
                end

                if nunk(i_line) -  c.h*FWHM_voigt(i_line) * n_FWHM_hi_res < nu(1) || nunk(i_line) +  c.h*FWHM_voigt(i_line) * n_FWHM_hi_res > nu(end)
                    continue
                end

                %% Voigt Binning Function
                switch binning_function
                    case 'SampleAtBinEdge'
                        % Sample only in the range available in the
                        % Holtsmark table H(u) (denoted edge)

                        % Sample high resolution only a few FWHM's from the
                        % center (denoted bnd). Otherwise low resolution

                        [u_low_res_mapped,nu_low_res_edge,nu_hi_res_sample,u_hi_res_mapped,len_phi,indcs]=split_low_hi_res_map_to_u(nu,nu_hi_res,nunk(i_line),Dnu(i_line),c.h*FWHM_voigt(i_line) * n_FWHM_hi_res,min_u,max_u,length(H_table.u),low_res_spacing);

                        if isempty(u_low_res_mapped)
                            continue
                        end

                        phi = zeros(1,len_phi);
                        phi(indcs.low_res_bin) = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_low_res_mapped,Mode);
                        phi_hi_res = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_hi_res_mapped,Mode);

                        % So as not too add phi twice in the hi-res region,
                        % we subtract a linear interpolation of low_res phi
                        % in the high res boundary
                        phi_hi_res = phi_hi_res - (phi(indcs.nu_low_bnd_edge) + ( nu_hi_res_sample - nu(indcs.nu_low_bnd) ) .* ( phi(indcs.nu_hi_bnd_edge) - phi(indcs.nu_low_bnd_edge) ) ./ ( nu(indcs.nu_hi_bnd) - nu(indcs.nu_low_bnd) ));
                        % Likewise, make the interpolant true for
                        % phi_low_res in this range
                        if indcs.nu_hi_bnd > indcs.nu_low_bnd +1
                            phi(indcs.nu_low_bnd_edge + 1 : indcs.nu_hi_bnd_edge - 1) = phi(indcs.nu_low_bnd_edge) + ( nu(indcs.nu_low_bnd + 1 : indcs.nu_hi_bnd - 1) - nu(indcs.nu_low_bnd) ) .* ( phi(indcs.nu_hi_bnd_edge) - phi(indcs.nu_low_bnd_edge) ) ./ ( nu(indcs.nu_hi_bnd) - nu(indcs.nu_low_bnd) );
                        end
                    case 'fConserving' % Oscillator Strength Conserving
                        [u_low_res_mapped,nu_low_res_edge,nu_hi_res_sample,u_hi_res_mapped,len_phi,indcs,nu_mid_hi_res_sample,u_mid_hi_res,u_mid_hi_res_mapped]=split_low_hi_res_map_to_u_fconsv(nu,nu_hi_res,nu_mid_hi_res,nunk(i_line),Dnu(i_line),c.h*FWHM_voigt(i_line) * n_FWHM_hi_res,c.h*FWHM_voigt(i_line)*Mode.kappa.n_FWHM_fcons,min_u,max_u,length(H_table.u),low_res_spacing);

                        if isempty(u_low_res_mapped)
                            continue
                        end

                        phi = zeros(1,len_phi);
                        phi(indcs.low_res_bin) = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_low_res_mapped,Mode);
                        phi_cum = sign(u_mid_hi_res).*interp_wrapper(interp_wrapper(H_table.H_int',a_mapped(i_line),Mode),u_mid_hi_res_mapped,Mode)/pi^(1/2);
                        phi_f_cons = diff(phi_cum)./diff(nu_mid_hi_res_sample);
                        phi_edge_sample = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_hi_res_mapped,Mode);
                        % Likewise, make the interpolant true for
                        % phi_low_res in this range
                        if indcs.nu_hi_bnd > indcs.nu_low_bnd +1
                            phi(indcs.nu_low_bnd_edge + 1 : indcs.nu_hi_bnd_edge - 1) = phi(indcs.nu_low_bnd_edge) + ( nu(indcs.nu_low_bnd + 1 : indcs.nu_hi_bnd - 1) - nu(indcs.nu_low_bnd) ) .* ( phi(indcs.nu_hi_bnd_edge) - phi(indcs.nu_low_bnd_edge) ) ./ ( nu(indcs.nu_hi_bnd) - nu(indcs.nu_low_bnd) );
                        end

                        phi_hi_res = zeros(size(indcs.hi_res_bin));
                        phi_hi_res(indcs.f_cons) =  phi_f_cons;
                        phi_hi_res(~indcs.f_cons) =  phi_edge_sample;
                        phi_hi_res = phi_hi_res - (phi(indcs.nu_low_bnd_edge) + ( nu_hi_res_sample - nu(indcs.nu_low_bnd) ) .* ( phi(indcs.nu_hi_bnd_edge) - phi(indcs.nu_low_bnd_edge) ) ./ ( nu(indcs.nu_hi_bnd) - nu(indcs.nu_low_bnd) ));
                    case 'DeltaOnly'
                        phi = 1./diff_nu(nunk_ind(i_line));
                end

                %% Add Voigt (phi) into kappa
                if length(phi)>1
                    %                     phi_norm_check(i_line)=trapz(nu,phi);
                    if Mode.kappa.add_omega4 % After Seaton 1994
                        phi = phi .* (( nu_low_res_edge./ nunk(i_line) ).^4.*(nu_low_res_edge<nunk(i_line))+(nu_low_res_edge>nunk(i_line)));
                        phi_hi_res = phi_hi_res .* (( nu_hi_res_sample ./ nunk(i_line) ).^4.*(nu_hi_res_sample<nunk(i_line))+(nu_hi_res_sample>nunk(i_line)));
                    end

                    dkappa_dsig = ij_no_frac*n_n(i_line).*(1-exp(-nunk(i_line)./T'))*w_hm_DAM_bb(i_line)/c.amu;
                    kappa_line = dkappa_dsig*signk(i_line);
                    if expansion_limit
                        kappa_line_hi_res = min(kappa_l_sob_lim.*nunk(i_line)/c.h,dkappa_dsig*signk(i_line));
                    else
                        kappa_line_hi_res = dkappa_dsig*signk(i_line);
                    end
                    k_bb(indcs.nu_low_edge:indcs.nu_hi_edge) = k_bb(indcs.nu_low_edge:indcs.nu_hi_edge) + kappa_line*phi*c.h;
                    k_bb_hi_res(indcs.hi_res_bin) = k_bb_hi_res(indcs.hi_res_bin) + kappa_line_hi_res*phi_hi_res*c.h;
                else
                    dkappa_dsig = ij_no_frac*n_n(i_line).*(1-exp(-nunk(i_line)./T')).*w_hm_DAM_bb(i_line)/c.amu;
                    if expansion_limit
                        kappa_line = min(kappa_l_sob_lim.*nunk(i_line)/c.h , signk(i_line).*dkappa_dsig);
                    else
                        kappa_line = signk(i_line).*dkappa_dsig;
                    end
                    k_bb_hi_res(nunk_ind(i_line)) = k_bb_hi_res(nunk_ind(i_line)) + kappa_line.*phi*c.h;
                end
            end


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
ExRateCoef = 0; %Will add this later...

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