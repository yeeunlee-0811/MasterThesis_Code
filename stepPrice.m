function [dhp] = stepPrice(inp)
%%%%input
XSinter = inp.XSinter;
price = XSinter(:,1);
eps = inp.eps;
J = inp.J;

% Computation of steps for the derivative of delta
dhp = zeros(J,1); % citydum : the number of cities : --> 1 in my case? // 3 --> 63?
for i =1:J
x0 = price(i); % 1 by 1
ax0 = abs(x0); % 1 by 1
if x0 ~= 0 
    dax0 = x0./ax0 ;
else 
    dax0 = 1;
end

dh1  = eps*max([ax0,(1e-2)*ones(size(x0,1),1)],[],2).*dax0;

xdh =x0+dh1;
dh1 = 1000*(xdh-x0);

dhp(i) =dh1;
end 

