function [Gam] = Gamma(inp, Dlikelmat)
% function Gamma calculates Gamma, the matrix of moment gradients w.r.t
% theta and delta
% Gamma : parmnum1 by parmnum1 (J + parmnum)

%%%%%%% Input 
%%% From INP
J =inp.J;
R = inp.R;

%%% From INPR
deltafin=inp.delta;
share0 = inp.S;
theta0 = inp.theta;
dhd = inp.dhd;
dh = inp.dh;
ntheta = inp.ntheta;
parmnum1 = inp.parmnum1;

Gam = zeros(parmnum1, parmnum1);
%%%%%%% Derivative of share moments
%% Derivative w.r.t delta

smldelta = deltafin;
theta = theta0;
ind = 1;
while ind <= J 
    chg = dhd(ind);
    smldelta(ind) = smldelta(ind) + chg;
    [mshare] =  predicted_share(theta,smldelta,inp); % J by 1
    
   fprintf('Share : Current number for delta deriv is %d\n',ind);
   Gam(1:J,ind) = (share0 - mshare).*(1/chg); 
   smldelta(ind) =smldelta(ind)-chg;
   ind = ind + 1;
end

%% Derivatiave w.r.t theta

theta= theta0;
k = 1;
while k <=ntheta
    chg = dh(k);
    theta(k) = theta(k)+chg;
    smldelta = deltafin;
    mshare = predicted_share(theta,smldelta,inp);
    
    fprintf('Share : Current number for theta deriv is %d\n',k);
    Gam(1:J,J+k) = (share0 - mshare).*(1/chg);
    theta(k)= theta(k)-chg;
    k = k+1;
end

%%%%%% Derivative of likelihood function moments
%% Derivative w.r.t delta

delta = deltafin ;
theta = theta0;
ind = 1;
while ind<=J
    chg = dhd(ind) ;
    delta(ind) = delta(ind) + chg ; % only delta has changed
     fprintf('lkl : Current number for delta deriv is %d\n',ind);
    [ll0] = likelihd(theta,delta,inp); % real likelihood fn estimate (not the FOC of likelihood fn)
    
    k=1;
    while k <= ntheta
        DIO = sum(Dlikelmat(:,k))/R; % delta is unchanged , only theta has changed
        chg1 = dh(k);
        theta(k) = theta(k) + chg1;
        fprintf('lkl : Current number for theta deriv under delta is %d\n',k);
        [ll1] = likelihd(theta,delta,inp); % delta and theta is changed
        
        Gam(J+k,ind) =((ll1-ll0)/chg1-DIO)/chg; 
        theta(k) = theta(k)-chg1;
        k = k+1;
    end
delta(ind)=delta(ind)-chg;
ind = ind+1;
end
    
%% REST 

%LLO
delta=deltafin;
theta = theta0;
[LL0]=likelihd(theta,delta,inp) %likelihood function

%%
j=1;
%% 
while j <= ntheta
    k=1;
    while k<=j
        chg = dh(k)
        theta(k) = theta(k)+chg
        
        fprintf('lkl : Current number for theta deriv is %d\n',k);
        [LL1] = likelihd(theta,delta,inp)      % LL1 : lklfn estimate if theta(k) is slightly changed
        
        chg2 = dh(j)              
        theta(j) = theta(j) + chg2
        
        fprintf('lkl : Current number for theta deriv under given theta is %d\n',j);
        [LL3] = likelihd(theta,delta,inp)      % LL3 : lklfn estimate if theta(j)&theta(k) is slightly changed
        
        theta(k)=theta(k)-chg
        
        [LL2] = likelihd(theta,delta,inp)      % LL2 : lklfn estimate if theta(j) is slightly changed 
        
        Gam(J+k,J+j) = ((LL3-LL2)-(LL1-LL0))/(chg*chg2*R)
        Gam(J+j,J+k) = ((LL3-LL2)-(LL1-LL0))/(chg*chg2*R)
        
        theta(j)=theta(j)-chg2;
        
        k=k+1;
    end
    j=j+1;
end


        