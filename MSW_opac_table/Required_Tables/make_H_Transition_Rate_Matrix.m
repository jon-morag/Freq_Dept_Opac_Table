function [out]=make_H_Transition_Rate_Matrix()
% build transition rate and oscillator strength matrices for Hydrogen levels decay rates

c=set_consts;

Z = 1;
I = c.IH;
a0 = c.a0; % bohr radius

% decay rates from excited states, up to n=150;
nmax = 150;
% nmax = 10;
n = 1:nmax;
Ank = zeros(length(n));

fnk = zeros(length(n));
Bnk = zeros(length(n));
nunk = zeros(length(n));

gn = 2*n.^2;

nprime = factorial(1:nmax);

for n=2:nmax
    k=1:n-1;
    nu_nk=I*(k.^-2-n^-2);
       
    cc = n-k;
    cprime = nprime(cc);
    kprime = nprime(k);
    
    x = -4*n*k./(n-k).^2;
    x2 = 1./x;

%     FH1=hypergeom([-n+1,-k],1,x);
%     FH2=hypergeom([-k+1,-n],1,x);
%     fnk1 = 32/3*k.^-2.*abs(((n-k)./(n+k)).^(2*n+2*k).*(n.*k).^-2.*(k.^-2-n.^-2).^-3.*(FH1+FH2).*(FH1-FH2)./(n-k)); % oscillator strength

    FH1 = 0*k;
    FH2 = 0*k;
    FH21 = 0*k;
    FH32 = 0*k;
    for kk=k
        FH1(kk)=hypergeom([-kk,-n+1],1,x(kk));
        FH2(kk)=hypergeom([-kk+1,-n],1,x(kk));
        FH21(kk) = hypergeom([-kk,-kk+1],cc(kk)+1,x2(kk));
        theta = kk*cc(kk)/(2*kk+cc(kk));
        FH32(kk) = hypergeom([-kk,-kk,theta+1],[cc(kk)+1,theta],x2(kk));
    end
    
    % old calculation, Menzel & Pekeris '35
    %     fnk1 = 32/3*n.^4.*k.^2.*(n-k).^(2*n+2*k-4).*(1./(n+k)).^(2*n+2*k+3).*(FH1+FH2).*(FH1-FH2); % oscillator strength
         fnk1 = 32/3*n.^4.*k.^2.*((n-k)./(n+k)).^(2*n+2*k+3).*(n-k).^(-7).*(FH1+FH2).*(FH1-FH2); % oscillator strength
       
    
    % new calculation according to Menzel '69
    fnk2 = (2./cc).^(2-2*cc).*k./3./cprime(k).^2.*(nprime(n)./kprime(k)./k.^cc(k)).^2.* ...
               ((1+cc./k).^(2*k+2)./(1+cc(k)/2./k).^(4*k+2*cc(k)+3)).*FH21.*FH32;

    fnktemp = [fnk1; fnk2];
    mask1 = ~isnan(fnktemp(1,:)) & ~isinf(fnktemp(1,:)) & ~(fnktemp(1,:)==0);
    mask2 = ~isnan(fnktemp(2,:)) & ~isinf(fnktemp(2,:)) & ~(fnktemp(2,:)==0);
    fnk(mask1,n) = fnktemp(1,mask1);
    fnk(mask2 & ~mask1,n) = fnktemp(2,mask2 & ~mask1);
    
    Ank(k,n) = gn(k)./gn(n)*2*c.c/a0.*(nu_nk/(c.me*c.c^2)).^2.*fnk(k,n)'; % spontaneous transition rate (Einstein A coefficient)
%%     Ank(n,k)=gn(k)./gn(n)*2*c.c/a0.*(nu_nk/(c.me*c.c^2)).^2.*fnk(n,k); % spontaneous transition rate (Einstein A coefficient)
%%     Bnk(k,n)=Ank(n,k).*gn(n)./gn(k)*(c.c^3*c.h^2)./(8*pi*nu_nk.^3);

    Bnk(k,n) = (4*pi)*fnk(k,n)'*pi*c.q^2/c.me/c.c./nu_nk;
    nunk(k,n) = nu_nk;
end
DecayRate = zeros(nmax,nmax);
DecayRate_up_only = zeros(1,nmax);
for n = 2:nmax
    DecayRate_up_only(n) = sum(Ank(:,n));
    for k = 1:n-1
        DecayRate(k,n) = sum(Ank(:,n))+sum(Ank(:,k));
%         DecayRate(k,n) = sum(Ank(n,:))+sum(Ank(k,:));
    end
end

out.Ank = Ank;
out.fnk = fnk;
out.Bnk = Bnk;
out.nunk = nunk;
out.DecayRate = DecayRate;
out.DecayRate_up_only = DecayRate_up_only;