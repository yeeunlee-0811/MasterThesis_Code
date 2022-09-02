function F  = find_mcallot(mc,inp,ownership,mctot)

% given mc, calculates optimal price
[optp,P_opt] = findp_cont(mc,inp,ownership);

P1 = P_opt(1:3)./sum(P_opt(1:3));
P2 = P_opt(4:6)./sum(P_opt(4:6));
P3 = P_opt(7:9)./sum(P_opt(7:9));
mc1 = mc(1:3);
mc2 = mc(4:6);
mc3 = mc(7:9);

F(1) = P1'*mc1 - mctot(1);
F(2) = P2'*mc2 - mctot(2);
F(3) = P3'*mc3 - mctot(3);
end