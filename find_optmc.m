function [F,mcspf,P_opt,pi,cpperson]=find_optmc(unitmc,discrate,ownership,inp)
% This function caculates optimal unit marginal cost of data per firm
% discount rate is given, and assume that having the same discount rates
% is optimal for each firm. 
R = inp.R;
totmc = inp.totmc;

% discrate is 1X2 vector
disc = [0,discrate];
discv = repmat(disc,1,3);
discvt = discv';

% unitmc is 1X3 vector : perfirm mc
mc_unit = zeros(1,9);
for i = 1:3
    mc_unit(1,3*(i-1)+1:3*(i-1)+3) = repmat(unitmc(i),1,3);
end
mc_unitt=mc_unit';

data = inp.data;
discvec = 1-discvt;

mcspf=mc_unitt.*discvec.*data;

[optp,P_opt] = findp_cont(mcspf,inp,ownership);

pFp = zeros(9,1);
pFp(1:3) = P_opt(1:3)./sum(P_opt(1:3));
pFp(4:6) = P_opt(4:6)./sum(P_opt(4:6));
pFp(7:9) = P_opt(7:9)./sum(P_opt(7:9));

cpperson = zeros(3,1);
cpperson(1) = totmc(1)/(sum(P_opt(1:3))*R);
cpperson(2) = totmc(2)/(sum(P_opt(4:6))*R);
cpperson(3) = totmc(3)/(sum(P_opt(7:9))*R);

F(1) = sum(pFp(1:3).*mcspf(1:3)) - cpperson(1);
F(2) = sum(pFp(4:6).*mcspf(4:6)) - cpperson(2);
F(3) = sum(pFp(7:9).*mcspf(7:9)) - cpperson(3);

pi = zeros(3,1);
pi(1) = sum((P_opt(1:3)*R).*(optp(1:3)-mcspf(1:3)));
pi(2) = sum((P_opt(4:6)*R).*(optp(4:6)-mcspf(4:6)));
pi(3) = sum((P_opt(7:9)*R).*(optp(7:9)-mcspf(7:9)));
pi= pi';
mcspf = mcspf';

end

