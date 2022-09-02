function [lklfn2,P,PS1,delta] = lklfn2_inter(inp,theta)

R=inp.R;
J=inp.J  ;
JN = inp.JN;
ninter1 = inp.ninter1;
ninter2 = inp.ninter2;
logS = inp.logS;

XIinter1=inp.XIinter1;
XIinter1=XIinter1';
XSinter1 = inp.XSinter1;

XIinter2=inp.XIinter2;
XIinter2=XIinter2';
XSinter2 = inp.XSinter2;

nup = inp.nup;
nudata =inp.nudata;
ssp = inp.ssp;
data = inp.data;

% allocating parameters
theta2 = theta;
if inp.inter1
    mu = 0;
    for k=1:ninter1 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter1(k,:).*repmat(XSinter1(:,k),1,R);
    end
    
else; inp.inter2;
    
    mu = 0;
    for k=1:ninter2 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter2(k,:).*repmat(XSinter2(:,k),1,R);
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