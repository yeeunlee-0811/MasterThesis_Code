function [delta,P] = cal_deltaP(theta,inp)
R=inp.R;
nikrc = inp.nikrc ;
J=inp.J  ;
JN = inp.JN;
ninter = inp.ninter;
logS = inp.logS;
nskrc = inp.nskrc;
XI=inp.XI ;
XS=inp.XS ;
XIinter=inp.XIinter;
XIinter=XIinter'; %XIinter changed to (ninter by R)
XSinter = inp.XSinter;
Xdp = inp.Xdp;


mu = 0;
for k=1:ninter % loop through all product characteristics except price
    mu = mu + theta(k)*XIinter(k,:).*repmat(XSinter(:,k),1,R);
end


tol = 1e-14;
dst = 1;
delta = 0.5*ones(N,1); % arbitrarily chosen starting values for delta
while dst > tol
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = sum(eU,1);
    P = mean(eU./denom,2);
    delta_new = delta + logS - log(P);
    dst = max(abs(delta - delta_new));
    delta = delta_new;
end

%if nargout>1 % if more than one output argument
    % calculate choice probabilities based on the final (or given) delta
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = sum(eU,1);
    P = mean(eU./denom,2);
end
