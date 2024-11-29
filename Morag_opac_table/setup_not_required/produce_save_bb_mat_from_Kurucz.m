function [Transition_Rate_Matrix] = produce_save_bb_mat_from_Kurucz(bb_lines,max_element,matrix_save_name)
c = set_consts();

ionization_state = round(mod(bb_lines.elem_chrg,1)*100);
element = floor(bb_lines.elem_chrg);
J1 = bb_lines.J1;
J2 = bb_lines.J2;

gfnk = 10.^(bb_lines.log_GF);

DecayRate = 10.^(bb_lines.log_Gam_rad);
DecayRate(isnan(DecayRate) | DecayRate == 1) = 0;
exratecoef = 10.^(bb_lines.log_Gam_Strk);
exratecoef(bb_lines.log_Gam_Strk == 0) = 0;
nu_nk = c.h*c.c ./ bb_lines.wl_nm * 1e7;

E1 = bb_lines.E1_cm_1*c.c*c.h;
E2 = bb_lines.E2_cm_1*c.c*c.h;
[Elo Elo_ind] = min([E1 E2],[],2);
Jlo = J1.*( Elo_ind == 1 ) + J2.*( Elo_ind == 2 );
Jhi = J2.*( Elo_ind == 1 ) + J1.*( Elo_ind == 2 );
glo = 2*Jlo+1;
ghi = 2*Jhi+1;
Ehi = E2.*( Elo_ind == 1 ) + E1.*( Elo_ind == 2 );
del_lambda = 1 - nu_nk./(Ehi-Elo);
tolerance = 0.03;
is_good = ~isnan(del_lambda) & abs(del_lambda) < tolerance;
is_bad = ~isnan(del_lambda) & abs(del_lambda) > tolerance;

%% Now redistribute these in a usable matrix
for ix = 1:min(max(element),max_element)
    for i_j = 0:ix-1
        ind = element == ix & ionization_state == i_j & ~is_bad; %Doing is_good keeps out values without an upper energy level.
        Transition_Rate_Matrix(ix,i_j+1).nu_nk = nu_nk(ind);
        Transition_Rate_Matrix(ix,i_j+1).element = element(ind);
        Transition_Rate_Matrix(ix,i_j+1).Ionization_state = ionization_state(ind);
        Transition_Rate_Matrix(ix,i_j+1).Elo = Elo(ind);
        Transition_Rate_Matrix(ix,i_j+1).Ehi = Ehi(ind);
        Transition_Rate_Matrix(ix,i_j+1).ExRateCoef = exratecoef(ind);
        Transition_Rate_Matrix(ix,i_j+1).DecayRate = DecayRate(ind);
        Transition_Rate_Matrix(ix,i_j+1).glo = glo(ind);
        Transition_Rate_Matrix(ix,i_j+1).ghi = ghi(ind);
        Transition_Rate_Matrix(ix,i_j+1).del_lambda = del_lambda(ind);
        Transition_Rate_Matrix(ix,i_j+1).fnk = gfnk(ind)./glo(ind);
    end
end

save(matrix_save_name,'Transition_Rate_Matrix','-v7.3')
end