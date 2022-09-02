function [res, theta3_0,theta1] = res(inp)
J=inp.J;
R=inp.R;
% the number of product characteristics
Krc=inp.Krc;
K=inp.K ;
XI=inp.XI ;
X=inp.X ;
JN=inp.JN ;
logS=inp.logS ;
nikrc=inp.nikrc ;
Z=inp.Z ;



%% residual 
Y= delta;
W = (Z'*Z)/J; % initial weighting matrix
xzWz = X'*Z*(W\Z');
    A = xzWz*X;
    B = xzWz*Y;
    theta1 = A\B;

res = Y - X*theta1;


end