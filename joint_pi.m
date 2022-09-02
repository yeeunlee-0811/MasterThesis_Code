function [jointp]=joint_pi(optp,inp,ownership)

ssp=optp;
Own = inp.Own;
cprice0 = optp; 
deltac = zeros(J,1);
alphac = theta1_indiv(2,:);
res2 = inp.res2;

theta2 = thetafin_indiv(1:ninter);
sigp = thetafin_indiv(ninter+1);
sigd = thetafin_indiv(ninter+2);

Exp1 = data; % drop original price
Exp122 = [ssp,Exp1];
Exp133 = [ones(J,1),Exp122]; % total 7 variables

delta_pchg = Exp133*theta1_indiv+res2; % Jx1

XSinter1 = repmat(cprice0,1,ninterp); % Jxninter
XSinter = [XSinter1,dt];

mu = 0;
for k=1:ninter % loop through all product characteristics except price
    mu = mu + theta2(k)*XIinter(k,:).*repmat(XSinter(:,k),1,R);
end
mu = mu + sigp*nup.*repmat(ssp,1,R);
mu = mu + sigd*nudata.*repmat(data,1,R);

U = repmat(delta_pchg,1,R) + mu;
eU = exp(U);
denom = 1+ sum(eU,1);
P = mean(eU./denom,2);
PS = eU./denom;

% J x 1 x R
PP = reshape(PS,J,1,R);
% J x J x R
PP = repmat(PP,1,J,1); % Create Nt identical copies of P and stack horizontally
PPt = permute(PP,[2 1 3]);
I = repmat(eye(J),1,1,R);

muc = 0;
for k=1:ninterp % loop through all product characteristics except price
    muc = muc + theta2(k)*XIinter(k,:);
end
muc = muc +sigp.*nup;

indiv = muc;
indiv = reshape(indiv,1,1,R);
indiv = repmat(indiv,J,1,1);
indiv = repmat(indiv,1,J,1);
alpha = repmat(alphac,J,J,R);

coeff = indiv +alpha;

dPc = mean(coeff.*PP.*(I - PPt),3); 

Ownmono = ones(9,9);
mc_org_own = cprice0 + inv((Own.*dPc)')*P;
% calculates joint profit
mc = optp + inv((Own.*dPc)')*P;
mc_org_mono = cprice0 +inv((Ownmono.*dPc)')*P;
pi = P_opt.*(optp-mc);
jointp=-sum(pi);
end
