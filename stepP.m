function [dh] = stepP(inp)

%%%%% INPUT
if inp.ponly
    theta = inp.theta1_p;
elseif inp.inter1
    theta = inp.theta1_it1;
elseif inp.wd2
    theta = inp.theta1_wd2;
elseif inp.wd3
    theta = inp.theta1_wd3;
end

eps = inp.eps;


x0=theta;
ax0 = abs(x0);
if x0 ~= 0
    dax0 = x0./ax0; % parmnum by 1
else 
    dax0 = 1;
end

dh = eps*max([ax0,(1e-2)*ones(size(x0,1),1)],[],2).*dax0;

xdh =x0 + dh;
dh =1000*(xdh-x0); % parmnum by 1
end