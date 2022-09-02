function piskt =  pi_skt(mcskt,mckt,mclg,discrate,ownership,inp)

% discrate is 1X2 vector
disc = [0,discrate];
discv = repmat(disc,1,3);
discvt = discv';

unitmc = [mcskt,mckt,mclg];
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
piskt = -sum(P_opt(1:3).*(optp(1:3)-mcspf(1:3)));