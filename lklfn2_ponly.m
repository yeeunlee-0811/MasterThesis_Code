function [lklfn2,P,PS1,delta] = lklfn2_ponly(inp,theta)

R=inp.R;
J=inp.J  ;
JN = inp.JN;
logS = inp.logS;
nup = inp.nup;
ssp = inp.ssp;
nudata = inp.nudata;
data = inp.data;
% allocating parameters
theta1= theta(1);
theta2 = theta(2);

mu = 0;
mu = mu + theta1*nup.*repmat(ssp,1,R);
mu = mu + theta2*nudata.*repmat(data,1,R);

tol = 1.0001e-10;
dst = 1;
delta = 0.5*ones(J,1); % arbitrarily chosen starting values for delta

while dst > tol
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = 1+ sum(eU,1);
    PS1=eU./denom;
    P = mean(eU./denom,2);
    %sum(isnan(P),'all');
    delta_new = delta + logS - log(P);
    % delta_new(1) = 0;
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