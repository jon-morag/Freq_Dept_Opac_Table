function H_table=HoltsmarkTable()
% Cumulative Holtsmark distribution
% the function ouputs the integral over beta = F/F0 for the Holtsmark distribution

Hfunc = @(beta) 2*beta./pi.*quadgk(@(y) y.*exp(-y.^(3/2)).*sin(beta.*y),0,inf);

beta = logspace(-6,log10(580),1000);

HH = 0*beta;
for ibeta = 1:length(beta)
    HH(ibeta) = Hfunc(beta(ibeta));
end

HH(beta>100) = 1.5*beta(beta>100).^(-5/2);

H = cumtrapz(beta,HH);

H_table.beta = beta;
H_table.H = H;

H_table.F = griddedInterpolant(beta,H);