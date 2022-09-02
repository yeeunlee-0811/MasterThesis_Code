function [res,theta1] = residfun(inp,theta)
% This function calculates residual which is used for GMM estimation

R =inp.R;
JN = inp.JN;
ninter = inp.ninter;
logS = inp.logS;
nskrc = inp.nskrc;
XS=inp.XS ;

Z=inp.zD;
W = inp.W;


if inp.conc % "concentrate out" linear parameters
    theta2 = theta;  
else
    theta1 = theta(1:nskrc);
    theta2 = theta(nskrc+1:end);
end

[~,~,~,delta] = lklfn2(inp,theta2);

Y = delta;

if inp.conc
    xzWz = XS'*Z*(W\Z');
    A = xzWz*XS;
    B = xzWz*Y;
    theta1 = A\B;
end

res = Y - XS*theta1;
end
