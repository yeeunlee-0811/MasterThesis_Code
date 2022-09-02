function [lklfn2,P,PS1,delta] = lklfn2_indiv(inp,theta)
% find estimates when there are individual terms 
% to do this divide theta

R=inp.R;
J=inp.J  ;
JN = inp.JN;
ninter = inp.ninter;
ninterp = inp.ninterp;
logS = inp.logS;
XSinter = inp.XSinter;
XIinter = inp.XIinter;
XIinter = XIinter';
nindiv = inp.nindiv;


mu = 0;
for k=1:ninter % loop through all product characteristics except price
    mu = mu + theta(k)*XIinter(k,:).*repmat(XSinter(:,k),1,R); % J by Rindiv
end

tol = 1.000050004108521e-10;
dst = 1;
delta = 0.5*ones(J,1); % arbitrarily chosen starting values for delta
delta(1) = 0;
while dst > tol
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    %sum(isnan(eU),'all');
    denom = sum(eU,1);
    %sum(isnan(denom),'all');
    PS1=eU./denom;
    P = mean(eU./denom,2);
    %sum(isnan(P),'all');
    delta_new = delta + logS - log(P);
    delta_new(1) = 0;
    dst = max(abs(delta - delta_new));
    delta = delta_new;
end

%if nargout>1 % if more than one output argument
    % calculate choice probabilities based on the final (or given) delta
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = sum(eU,1);
    P = mean(eU./denom,2);

PS1 = eU./denom ; 


logPS1 = log(PS1);
A = logPS1.*JN;

lklfn2 = -sum(sum(A));
end