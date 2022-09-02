function [lklfn2,P,PS1,delta] = lklfn2_blp(inp,theta)

R=inp.R;
J=inp.J  ;
JN = inp.JN;
ninter = inp.ninter;
logS = inp.logS;
nskrc = inp.nskrc;
XIinter=inp.XIinter;
XIinter=XIinter';
XSinter = inp.XSinter;
nup = inp.nup;
nudata =inp.nudata;
ssp = inp.ssp;
data = inp.data;
income = inp.income;


% allocating parameters
thetap = theta;

mu = 0;
mu = mu + (thetap./income).*repmat(ssp,1,R);


tol = 1.0001e-10;
dst = 1;
delta = 0.5*ones(J,1); % arbitrarily chosen starting values for delta
% delta(1) = 0;
while dst > tol
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    %sum(isnan(eU),'all');
    denom = 1+ sum(eU,1);
    %sum(isnan(denom),'all');
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