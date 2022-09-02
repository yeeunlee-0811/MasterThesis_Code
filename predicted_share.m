function [P,PS,U] = predicted_share(theta,delta,inp)
% This function calculates predicted products shares and individual choice
% probabilities

R=inp.R;
ninter1 = inp.ninter1;
ninter2 = inp.ninter2;
ninter3 = inp.ninter3;
ninter4 = inp.ninter4;
ninter5 = inp.ninter5;

XIinter1=inp.XIinter1;
XIinter1 =XIinter1';
XIinter2 = inp.XIinter2;
XIinter2 = XIinter2';
XIinter3 = inp.XIinter3;
XIinter3 = XIinter3';
XIinter4 = inp.XIinter4;
XIinter4 = XIinter4';


XIinter5=inp.XIinter5;
XIinter5=XIinter5';
XSinter5 = inp.XSinter5;

XIinter6=inp.XIinter6;
XIinter6=XIinter6';
XSinter6 = inp.XSinter6;


XSinter1 = inp.XSinter1;
XSinter2 = inp.XSinter2;
XSinter3 = inp.XSinter3;
XSinter4 = inp.XSinter4;

nup = inp.nup;
nudata =inp.nudata;
ssp = inp.ssp;
data = inp.data;

% allocating parameters
if inp.ponly
    thetap = theta(1);
    thetad = theta(2);
    mu = 0;
    mu = mu + thetap*nup.*repmat(ssp,1,R);
    mu = mu + thetad*nudata.*repmat(data,1,R);
elseif inp.inter1
    theta1 =theta;
    mu = 0;
    for k=1:ninter1 % loop through all product characteristics except price
        mu = mu + theta1(k)*XIinter1(k,:).*repmat(XSinter1(:,k),1,R);
    end
elseif inp.wd2
    theta2=theta;
    mu = 0;
    for k=1:ninter4 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter4(k,:).*repmat(XSinter4(:,k),1,R);
    end
    elseif inp.wd3
    theta2=theta;
    mu = 0;
    for k=1:ninter5 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter5(k,:).*repmat(XSinter5(:,k),1,R);
    end
    elseif inp.wd4
    theta2=theta;
    mu = 0;
    for k=1:ninter5 % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter6(k,:).*repmat(XSinter6(:,k),1,R);
    end
end
%% 
U = repmat(delta,1,R) + mu;
eU = exp(U);
denom = 1+ sum(eU,1);
P = mean(eU./denom,2);
PS = eU./denom;


end