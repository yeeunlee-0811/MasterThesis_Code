function [lklfn2,P,PS1,delta] = lklfn2_wd(inp,theta)

R=inp.R;
J=inp.J  ;
JN = inp.JN;
ninter3 = inp.ninter3;
ninter4 = inp.ninter4;
ninter5 = inp.ninter5;

ninterp3 = inp.ninterp3;
ninterp4 = inp.ninterp4;
logS = inp.logS;

XIinter3=inp.XIinter3;
XIinter3=XIinter3';
XSinter3 = inp.XSinter3;

XIinter4=inp.XIinter4;
XIinter4=XIinter4';
XSinter4 = inp.XSinter4;

XIinter5=inp.XIinter5;
XIinter5=XIinter5';
XSinter5 = inp.XSinter5;

XIinter6=inp.XIinter6;
XIinter6=XIinter6';
XSinter6 = inp.XSinter6;

nup = inp.nup;
nudata =inp.nudata;
ssp = inp.ssp;
data = inp.data;

% allocating parameters
theta2 = theta;
if inp.wd2
    
    mu = 0;
    for k=1:ninter4 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter4(k,:).*repmat(XSinter4(:,k),1,R);
    end
elseif inp.wd3
    mu = 0;
    for k=1:ninter5 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter5(k,:).*repmat(XSinter5(:,k),1,R);
    end
    elseif inp.wd4
    mu = 0;
    for k=1:ninter5 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter6(k,:).*repmat(XSinter6(:,k),1,R);
    end
end

tol = 1.0001e-10;
dst = 1;
delta = 0.5*ones(J,1); % arbitrarily chosen starting values for delta

while dst > tol
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = 1+ sum(eU,1);
    PS1=eU./denom;
    P = mean(eU./denom,2);
    delta_new = delta + logS - log(P);
    dst = max(abs(delta - delta_new));
    delta = delta_new;
end

U = repmat(delta,1,R) + mu;
eU = exp(U);
denom = 1+ sum(eU,1);
P = mean(eU./denom,2);

PS1 = eU./denom ;


logPS1 = log(PS1);
A = logPS1.*JN;

lklfn2 = -sum(sum(A));
end