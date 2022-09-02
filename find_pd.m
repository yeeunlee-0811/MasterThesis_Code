function [fn_pdd]=find_pd(pd,inp)
data = inp.data;
S = inp.S;
for k = 1:3
    pdd((k-1)*3+1:3*k) = pd(k);
end
fn_pdd=pdd*data;
end