function [V1, Dlikelmat0, Hessian]=varMOM(inp)
%%%%% INPUT
S = inp.S;
R =inp.R;

delta0 = inp.delta;
theta0 = inp.theta;
ntheta=inp.ntheta;
parmnum1 = inp.parmnum1;
Dlikelmat0 = inp.Dlikelmat0; % Dlikelmat: obs by parmnum
dh=inp.dh;

%% Likelihood moment 
k=1;
while k <= ntheta 
    Blikel = dlkl(delta0,theta0,inp,k); % the sec moment cond. evaluated at optimum theta1 : R by 1
    chg = dh(k);
    theta0(k)=theta0(k)+chg;
    fprintf('Current number is %d\n',k);
    Dlikelmat0(:,k) = (dlkl(delta0,theta0,inp,k)- Blikel)./chg
    theta0(k) = theta0(k)-chg;
    k = k+1;
end 

%% Share moment
shrvec = S'; % observed share of products (1 by J)

%% V1
V1 = zeros(parmnum1,parmnum1) ; % J + parmnum = parmnum1
n = 1;
while n<=R
    shrind = [(-shrvec), Dlikelmat0(n,:)]; % 1 by J+parmnum
    V1 = V1 + shrind.*shrind' ;
    n = n+1 ;
end
V1=V1./R;
end
