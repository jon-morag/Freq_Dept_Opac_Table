function [k_bf_bb]=kappa_abs_bf_bb(T,rho,nu,Mode,Y,e_pop)
% calculates bf and bb for atomic mix
% bf: Analytic for H-like (one electron), Verner Yakovlev table otherwise
% bb: Analytic for H-like, Kurucz otherwise
% for efficiency, bb line peaks are computed in hi-res, and line wings are
% computed with lower freq. res 

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

%% bound-free opacity

natomm5n0 = rho*sum(X./A)/c.amu;
rion = (4*pi*natomm5n0/3).^(-1/3);
F0 = Y*sum(X.*Z./A)/sum(X./A)./rion.^2;

num_frac = X./A/sum(X./A);

if Mode.kappa.bf_on
    k_bf = kappa_abs_bf(n, T, nu, Mode,e_pop,Z,num_frac,F0,c);
else
    k_bf = 0;
end

    %% Bound - bound opacity
nem5n0 = Y.*rho*(Z./A*X')/c.amu;
H_table = Mode.HjertingTable;     

plot_i = 1;
if Mode.kappa.bb_on
    [k_bb, k_bb_hi_res] = kappa_abs_bb(nu_hi_res,nu,T,Mode,F0,e_pop, Z,A, X./A/sum(X./A));
else
    k_bb = 0;
    k_bb_hi_res = 0;
end

k_bb_hi_res(k_bb_hi_res<0) = 0; %Due to the subtraction we perform, there can be a negative value due to poor interp.
k_bf_bb = interp1( nu , k_bf + k_bb , nu_hi_res) + k_bb_hi_res;

end

function k_bf = kappa_abs_bf(n, T, nu, Mode,e_pop,Z,num_frac,F0,c)
k_bf = 0;
for ix =  Mode.kappa.bf_species_select % select specific species
    if Mode.kappa.plot_n_eff && Z(ix)>=2; figure(ix) ; end
    for ij = 1:Z(ix)
        if Z(ix) == ij && Mode.kappa.bf_Hydrogenic_on %H-like
            k_bf = k_bf + kappa_abs_bf_ij_Hlike(n, T, nu, Mode,e_pop(ix),Z(ix),num_frac(ix),F0,c);
        elseif Z(ix) > ij && Mode.kappa.bf_VY_on
            k_bf = k_bf + kappa_abs_bf_ij_VY(T,nu, Mode,e_pop(ix),Z(ix),num_frac(ix),F0,c,ij);
        end
    end
    if Mode.kappa.plot_n_eff && Z(ix)>=2
        ax = gca(); ax.FontSize = 14;
        xlabel('State Number'); ylabel('State Energy');
        title(['NIST listing vs. Effective "H-like" energy - Z = ' num2str(Z(ix))]);
        legend(p_En,cellfun(@(ij) ['Ionization ' num2str(ij)] , num2cell(2:Z(ix)),'UniformOutput',0))
    end
end
end

function k_bf_ix_Hlike = kappa_abs_bf_ij_Hlike(n, T, nu, Mode,e_pop_i , Z_i , num_frac_i,F0,c)
IH=c.IH;
gn_lo = 2*n.^2;

En = IH*Z_i^2./n.^2;
if Mode.kappa.include_w_hm_Hydrogenic
    w_hm = get_w_HM(F0,n,1,1,Mode,c,Z_i);
    Zr = gn_lo.*w_hm.*exp(-(IH*Z_i^2-En)./T');
else
    Zr = gn_lo.*exp(-(IH*Z_i^2-En)./T');
end
Gr = sum(Zr,2);
n_n = Zr./Gr;
nu_m_En = nu-En'; %+Icm5n0';
En_o_nu = En'./nu;
ind_allowed = (nu_m_En>0);

if Mode.kappa.Hydrogenic_Dappen_HM_factor
    % The equivalent n state
    n_DAM = real(( 1./n.^2 - nu/IH/Z_i.^2 ).^(-1/2));
    w_hm_DAM_bf = get_w_HM(F0,n_DAM,1,1,Mode,c,Z_i);

    nu_I = (nu_m_En).*(ind_allowed)+1.*(~ind_allowed);
    z = sqrt(En'./nu_I).*(ind_allowed)+1.*(~ind_allowed);
    g_nf = 8*pi*sqrt(3)*En_o_nu.*exp(-4*z.*acot(z))./(1-exp(-2*pi*z));
    sig_bf = n'./Z_i^2*8/sqrt(3)*c.ale^-3*c.sigT.*En_o_nu.^3;
    sig_bf = sig_bf .* (g_nf.*(nu_m_En>0)+8*pi*sqrt(3)*En_o_nu*exp(-4).*(nu_m_En<=0)); %This is a cheat, I use the gaunt factor value at zero
    sig_bf = sig_bf.*((1-w_hm_DAM_bf).*(nu_m_En<0) + 1.*(nu_m_En>=0));
else
    nu_I = (nu_m_En).*(ind_allowed)+1.*(~ind_allowed);
    z = sqrt(En'./nu_I).*(ind_allowed)+1.*(~ind_allowed);
    g_nf = 8*pi*sqrt(3)*En_o_nu.*exp(-4*z.*acot(z))./(1-exp(-2*pi*z));
    sig_bf = n'./Z_i^2*8/sqrt(3)*c.ale^-3*c.sigT.*En_o_nu.^3.*(g_nf.*(nu_m_En>0)+8*pi*sqrt(3)*En_o_nu*exp(-4).*(nu_m_En==0)+0.*(nu_m_En<0));
end

x = nu./T;
indlow = x<1e-8;

k_bf_ix_Hlike = e_pop_i.z(Z_i)*sum(n_n'.*sig_bf.*((1-exp(-x)).*(~indlow)+x.*(indlow)))*num_frac_i/c.amu;
end

function k_bf_ix_VY = kappa_abs_bf_ij_VY(T,nu, Mode,e_pop_i,Z_i,num_frac_i,F0,c,ij)
IH=c.IH;

[sigma_VY_bf_arr, n_VY] = sigma_VY_bf(nu,Z_i,ij,Mode.kappa.bf_VY_tbl,c);
if ~isempty(sigma_VY_bf_arr)
    In = Mode.Ionization.energies(Z_i,ij)*c.eV;
    gn = Mode.Ionization.levels(Z_i,ij).g; En=Mode.Ionization.levels(Z_i,ij).energies*c.eV;

    % Cut the levels by half to check if we are far from convergence
    levels_cutoff_length = ceil(length(gn)*Mode.kappa.levels_cutoff_frac);
    gn = gn(1:levels_cutoff_length);
    En = En(1:levels_cutoff_length);

    if Mode.kappa.include_w_hm_nonHydrogenic
        n_eff= real(ij*sqrt(IH./(In-En)));
        w_hm = get_w_HM(F0,n_eff,1,1,Mode,c,ij);
        Zr = gn.*w_hm.*exp(-En./T');
        Zr(find(n_eff==0)) = 0; %this means it was an imaginary result.
    else
        Zr = gn.*exp(-En./T');
    end

    Gr = sum(Zr,1);

    if Mode.kappa.extend_partition_func_n % && Zr(find(Zr,1,'last'))/Gr > Mode.kappa.e_partition_frac_cutoff %Turns out that this criterion is not always right. So always extend the partition function.
        if ~Mode.kappa.include_w_hm_nonHydrogenic; n_eff = ij * real(sqrt(IH./(In-En))); end % Don't calculate twice
        n_extend = round(max(n_eff)) + 1 : Mode.kappa.extension_quantum_number;
        %                         n_check = 1:round(max(n_eff));

        % For consistent addition assume the NIST ground state is zero.
        En_extend = In - ij^2*c.IH./n_extend.^2;

        if Mode.kappa.include_w_hm_nonHydrogenic
            w_hm_extend = get_w_HM(F0,n_extend,1,1,Mode,c,ij);
            Gr_extend = sum(2*n_extend.^2.*w_hm_extend.*exp(-En_extend./T'));
        else
            Gr_extend = sum(2*n_extend.^2.*exp(-En_extend./T'));
        end

        Gr = Gr + Gr_extend;
    end

    n_ground = sum(Zr(En<=0.3*c.eV))/Gr;
    %                     n_n = Zr./Gr;
    %All cross sections are given in terms of the ground
    %states.
    k_bf_ix_VY = e_pop_i.z(ij)*n_ground*sum(sigma_VY_bf_arr,1)*num_frac_i/c.amu;  

    if Mode.kappa.plot_n_eff
        % Examine the n_eff prescription
        p_En(ij) = plot(En/c.eV);
        hold on
        ax = gca(); ax.ColorOrderIndex = mod(ax.ColorOrderIndex-2,7)+1;
        plot(In*(1-1./round(n_eff).^2)/c.eV,'--')
        ax = gca(); ax.ColorOrderIndex = mod(ax.ColorOrderIndex-2,7)+1;
        len_n_eff = length(n_eff(n_eff>0));
        plot([len_n_eff len_n_eff+10] , [In/c.eV In/c.eV] ,'LineWidth',2.5)
    end
end

end

function [k_bb,k_bb_hi_res] = kappa_abs_bb(nu_hi_res,nu,T,Mode,F0,e_pop,Z,A,num_frac)
c=set_consts();
e_cutoff = 1e-4*c.eV; % Minimum transition energy

diff_nu = diff2(nu_hi_res);
nu_mid_hi_res = [nu_hi_res(1) , (nu_hi_res(2:end) + nu_hi_res(:,1:end-1,:))/2 , nu_hi_res(end)];
H_table = Mode.HjertingTable;
max_a = log(max(H_table.a));
min_a = log(min(H_table.a));
max_u = log(max(H_table.u));
min_u = log(min(H_table.u));
max_nu = log(max(nu));
min_nu = log(min(nu));

k_bb = zeros(1,length(nu));
k_bb_hi_res = zeros(1,length(nu_hi_res));

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

            if Z(ix) == ij && Mode.kappa.H_e_impact_broadening %H-like!
                %Full half - width %Armstrong 1964
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
            a_mapped = (length(H_table.a)-1)/(max_a-min_a)*(log(a)-min_a) + 1;

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
            nunk_ind = round((length(nu_hi_res)-1)/(max_nu-min_nu)*(log(nunk)-min_nu) + 1);

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
                        nu_low_edge_idx = find(nu > nunk(i_line)-max(H_table.u)/Dnu(i_line),1,'first');
                        nu_low_bnd_edge_idx = find(nu(nu_low_edge_idx:end) < (nunk(i_line) -  c.h*FWHM_voigt(i_line) * n_FWHM_hi_res),1,'last');
                        nu_low_bnd_idx = nu_low_edge_idx + nu_low_bnd_edge_idx-1;
                        nu_hi_bnd_idx = nu_low_bnd_idx + find(nu(nu_low_bnd_idx+1:end) > (nunk(i_line) +  c.h*FWHM_voigt(i_line) * n_FWHM_hi_res),1,'first');
                        if isempty(nu_low_bnd_idx) || isempty(nu_hi_bnd_idx)
                            continue
                        end
                        nu_hi_edge_idx = nu_hi_bnd_idx + find(nu(nu_hi_bnd_idx+1:end) < nunk(i_line)+max(H_table.u)/Dnu(i_line),1,'last');
                        nu_hi_bnd_edge_idx = nu_hi_bnd_idx - nu_low_edge_idx+1;
                        phi = zeros(1,1+nu_hi_edge_idx-nu_low_edge_idx);

                        low_res_bin_idx = [1:(1+nu_low_bnd_idx-nu_low_edge_idx) (1+nu_hi_bnd_idx-nu_low_edge_idx):length(phi)];
                        nu_low_res_edge = nu(nu_low_edge_idx:nu_hi_edge_idx);
                        nu_low_res_sample = nu([nu_low_edge_idx:nu_low_bnd_idx nu_hi_bnd_idx:nu_hi_edge_idx]);

                        nu_hi_res_low_idx = 2 + (nu_low_bnd_idx - 1)*low_res_spacing;
                        nu_hi_res_hi_idx = (nu_hi_bnd_idx - 1)*low_res_spacing;
                        hi_res_bin_idx = nu_hi_res_low_idx:nu_hi_res_hi_idx;
                        nu_hi_res_sample = nu_hi_res(hi_res_bin_idx);

                        u_low_res = (nu_low_res_sample-nunk(i_line)).*Dnu(i_line);
                        %u_low_res(abs(u_low_res)<min(H_table.u)) = min(H_table.u); % We don't care on which side of zero it is. Not sensitive.
                        u_low_res_mapped = (length(H_table.u)-1)/(max_u-min_u)*(log(abs(u_low_res))-min_u) + 1;

                        u_hi_res = (nu_hi_res_sample-nunk(i_line)).*Dnu(i_line);
                        u_hi_res_mapped = (length(H_table.u)-1)/(max_u-min_u)*(log(abs(u_hi_res))-min_u) + 1;
                        u_hi_res_mapped(u_hi_res_mapped>length(H_table.u)) = length(H_table.u);
                        u_hi_res_mapped(u_hi_res_mapped<1) = 1;

                        phi(low_res_bin_idx) = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_low_res_mapped,Mode);
                        phi_hi_res = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_hi_res_mapped,Mode);

                        % So as not too add phi twice in the hi-res region,
                        % we subtract a linear interpolation of low_res phi
                        % in the high res boundary
                        phi_hi_res = phi_hi_res - (phi(nu_low_bnd_edge_idx) + ( nu_hi_res_sample - nu(nu_low_bnd_idx) ) .* ( phi(nu_hi_bnd_edge_idx) - phi(nu_low_bnd_edge_idx) ) ./ ( nu(nu_hi_bnd_idx) - nu(nu_low_bnd_idx) ));
                        % Likewise, make the interpolant true for
                        % phi_low_res in this range
                        if nu_hi_bnd_idx > nu_low_bnd_idx +1
                            phi(nu_low_bnd_edge_idx + 1 : nu_hi_bnd_edge_idx - 1) = phi(nu_low_bnd_edge_idx) + ( nu(nu_low_bnd_idx + 1 : nu_hi_bnd_idx - 1) - nu(nu_low_bnd_idx) ) .* ( phi(nu_hi_bnd_edge_idx) - phi(nu_low_bnd_edge_idx) ) ./ ( nu(nu_hi_bnd_idx) - nu(nu_low_bnd_idx) );
                        end
                    case 'fConserving' % Oscillator Strength Conserving
                        nu_low_edge_idx = find(nu > nunk(i_line)-max(H_table.u)/Dnu(i_line),1,'first');
                        nu_low_bnd_edge_idx = find(nu(nu_low_edge_idx:end) < (nunk(i_line) -  c.h*FWHM_voigt(i_line) * n_FWHM_hi_res),1,'last');
                        nu_low_bnd_idx = nu_low_edge_idx + nu_low_bnd_edge_idx-1;
                        nu_hi_bnd_idx = nu_low_bnd_idx + find(nu(nu_low_bnd_idx+1:end) > (nunk(i_line) +  c.h*FWHM_voigt(i_line) * n_FWHM_hi_res),1,'first');
                        if isempty(nu_low_bnd_idx) || isempty(nu_hi_bnd_idx)
                            continue
                        end
                        nu_hi_edge_idx = nu_hi_bnd_idx + find(nu(nu_hi_bnd_idx+1:end) < nunk(i_line)+max(H_table.u)/Dnu(i_line),1,'last');
                        nu_hi_bnd_edge_idx = nu_hi_bnd_idx - nu_low_edge_idx+1;
                        phi = zeros(1,1+nu_hi_edge_idx-nu_low_edge_idx);

                        low_res_bin_idx = [1:(1+nu_low_bnd_idx-nu_low_edge_idx) (1+nu_hi_bnd_idx-nu_low_edge_idx):length(phi)];
                        nu_low_res_edge = nu(nu_low_edge_idx:nu_hi_edge_idx);
                        nu_low_res_sample = nu([nu_low_edge_idx:nu_low_bnd_idx nu_hi_bnd_idx:nu_hi_edge_idx]);

                        nu_hi_res_low_idx = 2 + (nu_low_bnd_idx - 1)*low_res_spacing;
                        nu_hi_res_hi_idx = (nu_hi_bnd_idx - 1)*low_res_spacing;
                        hi_res_bin_idx = nu_hi_res_low_idx : nu_hi_res_hi_idx;
                        hi_res_mid_bin_idx = (nu_hi_res_low_idx-1) : nu_hi_res_hi_idx;

                        nu_hi_res_sample = nu_hi_res( hi_res_bin_idx );
                        f_cons_idx = abs(nu_hi_res_sample-nunk(i_line)) <  c.h*FWHM_voigt(i_line)*Mode.kappa.n_FWHM_fcons;
                        nu_hi_res_sample_f_cons = nu_hi_res_sample(~f_cons_idx);
                        nu_mid_hi_res_sample = nu_mid_hi_res(hi_res_mid_bin_idx);
                        nu_mid_hi_res_sample = nu_mid_hi_res_sample(find(f_cons_idx + [0 f_cons_idx(1:end-1)]));
                        if f_cons_idx(end)
                            nu_mid_hi_res_sample = [nu_mid_hi_res_sample nu_mid_hi_res(nu_hi_res_hi_idx)];
                        end

                        u_low_res = (nu_low_res_sample-nunk(i_line)).*Dnu(i_line);
                        u_low_res(abs(u_low_res)<min(H_table.u)) = min(H_table.u); % We don't care on which side of zero it is. Not sensitive.
                        u_low_res_mapped = (length(H_table.u)-1)/(max_u-min_u)*(log(abs(u_low_res))-min_u) + 1;

                        u_hi_res = (nu_hi_res_sample_f_cons-nunk(i_line)).*Dnu(i_line);
                        u_hi_res_mapped = (length(H_table.u)-1)/(max_u-min_u)*(log(abs(u_hi_res))-min_u) + 1;
                        u_hi_res_mapped(u_hi_res_mapped<1) = 1;

                        u_mid_hi_res = (nu_mid_hi_res_sample-nunk(i_line)).*Dnu(i_line);
                        u_mid_hi_res_mapped = (length(H_table.u)-1)/(max_u-min_u)*(log(abs(u_mid_hi_res))-min_u) + 1;
                        u_mid_hi_res_mapped(u_mid_hi_res_mapped<1) = 1;

                        phi(low_res_bin_idx) = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_low_res_mapped,Mode);
                        phi_cum = sign(u_mid_hi_res).*interp_wrapper(interp_wrapper(H_table.H_int',a_mapped(i_line),Mode),u_mid_hi_res_mapped,Mode)/pi^(1/2);
                        phi_f_cons = diff(phi_cum)./diff(nu_mid_hi_res_sample);
                        phi_edge_sample = Dnu(i_line)/pi^(1/2).*interp_wrapper(interp_wrapper(H_table.H',a_mapped(i_line),Mode),u_hi_res_mapped,Mode);
                        % Likewise, make the interpolant true for
                        % phi_low_res in this range
                        if nu_hi_bnd_idx > nu_low_bnd_idx +1
                            phi(nu_low_bnd_edge_idx + 1 : nu_hi_bnd_edge_idx - 1) = phi(nu_low_bnd_edge_idx) + ( nu(nu_low_bnd_idx + 1 : nu_hi_bnd_idx - 1) - nu(nu_low_bnd_idx) ) .* ( phi(nu_hi_bnd_edge_idx) - phi(nu_low_bnd_edge_idx) ) ./ ( nu(nu_hi_bnd_idx) - nu(nu_low_bnd_idx) );
                        end

                        phi_hi_res = zeros(size(hi_res_bin_idx));
                        phi_hi_res(f_cons_idx) =  phi_f_cons;
                        phi_hi_res(~f_cons_idx) =  phi_edge_sample;
                        phi_hi_res = phi_hi_res - (phi(nu_low_bnd_edge_idx) + ( nu_hi_res_sample - nu(nu_low_bnd_idx) ) .* ( phi(nu_hi_bnd_edge_idx) - phi(nu_low_bnd_edge_idx) ) ./ ( nu(nu_hi_bnd_idx) - nu(nu_low_bnd_idx) ));
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

                    sig_bb = signk(i_line)*phi*c.h;
                    sig_bb_hi_res = signk(i_line)*phi_hi_res*c.h;
                    dkappa_dsig = ij_no_frac*n_n(i_line).*(1-exp(-nunk(i_line)./T'))*w_hm_DAM_bb(i_line)/c.amu;
                    k_bb(nu_low_edge_idx:nu_hi_edge_idx) = k_bb(nu_low_edge_idx:nu_hi_edge_idx) + dkappa_dsig*sig_bb;
                    k_bb_hi_res(hi_res_bin_idx) = k_bb_hi_res(hi_res_bin_idx) + dkappa_dsig*sig_bb_hi_res;
                else
                    sig_bb = signk(i_line)*phi*c.h;
                    k_bb_hi_res(nunk_ind(i_line)) = k_bb_hi_res(nunk_ind(i_line)) + ij_no_frac*n_n(i_line).*(1-exp(-nunk(i_line)./T')).*sig_bb.*w_hm_DAM_bb(i_line)/c.amu;
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

%% Artificially extend the partition function to check our sensitivity to finite energy levels
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