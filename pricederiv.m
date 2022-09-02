function [foc,dP,P] = pricederiv(price,inp)


% used for price change analysis 
%use estimated coefficients to calculate choice probabilites
theta1=inp.theta1;
theta2 = inp.theta;
res2= inp.res2; % Jx1
XSred = inp.XSred; % Explanatory variables without constant (6 variables)
XIinter = inp.XIinter;
XIinter= XIinter';
alpha = theta1(2,:);
J = inp.J;
ninter= inp.ninter;
R = inp.R;
Own=inp.Own;

Exp1 = XSred(:,2:end); % drop original price
Exp1 = [price,Exp1];
Exp1 = [ones(J,1),Exp1]; % total 7 variables

delta_pchg = Exp1*theta1+res2; % Jx1

XSinter = repmat(price,1,ninter); % Jxninter

mu = 0;
for k=1:ninter % loop through all product characteristics except price
    mu = mu + theta2(k)*XIinter(k,:).*repmat(XSinter(:,k),1,R);
end

U = repmat(delta_pchg,1,R) + mu;
eU = exp(U);
denom = sum(eU,1);
P = mean(eU./denom,2);
PS = eU./denom;

% J x 1 x R
PP = reshape(PS,J,1,R);
% J x J x R
PP = repmat(PP,1,J,1); % Create Nt identical copies of P and stack horizontally
PPt = permute(PP,[2 1 3]);
I = repmat(eye(J),1,1,R);

mu = 0;
for k=1:ninter % loop through all product characteristics except price
    mu = mu + theta2(k)*XIinter(k,:);
end
indiv=mu;
indiv = reshape(indiv,1,1,R);
indiv = repmat(indiv,J,1,1);
indiv = repmat(indiv,1,J,1);

coeff = indiv+alpha;

dP = mean(coeff.*PP.*(I - PPt),3); % NtxNt matrix of derivatives wrt price

foc = ((Own.*dP)')*(price) + P; % + penalty;
end

