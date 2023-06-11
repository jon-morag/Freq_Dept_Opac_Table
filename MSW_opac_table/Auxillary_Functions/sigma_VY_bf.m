% Analytic fits for partial photoionization cross sections.
% D.A. Verner and D.G. Yakovlev
function [sig_out,n_out] = sigma_VY_bf(nu,Z,ij,VY_tbl, c)
% %For He I:
% E_th = 24.59*c.eV; nu_0 = 5.996*c.eV; sig_0=4470*10^(-18); % Mb to cm^2
% y_a=2.1990; P = 6.098; y_W = 0; l = 0; 

VY_lines = VY_tbl(VY_tbl.NZ == Z & VY_tbl.Ion_State == Z-ij+1,:);

if ~isempty(VY_lines)
    for i_VY = 1:height(VY_lines)
        VY_line = VY_lines(i_VY,:);
        E_th = VY_line.EionizationThreshold*c.eV; nu_0 = VY_line.E_0*c.eV; l = VY_line.lorbitalqn;
        sig_0 = VY_line.sigma_0*1e-18; y_a = VY_line.y_a; P = VY_line.P; y_W = VY_line.y_W;
        
        y = nu/nu_0; Q = 5.5 + l - 0.5*P;
        sig = sig_0 *( (y-1).^2 + y_W.^2 ) .* y.^(-Q) .* (1+sqrt(y/y_a)).^(-P) ;
        sig_out(i_VY,:) = sig.*(nu >= E_th);
        n_out(i_VY,:) = VY_line.nquantumnumber;
    end
else
    sig_out = [];
end

end