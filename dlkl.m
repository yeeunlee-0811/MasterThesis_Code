function [dlklfn] = dlkl(delta0,theta0,inp,k)

%%%%% INPUT

JN=inp.JN;
R = inp.R;


step = 1e-5; % step-length for forward-backward finite differences
dlklfn = zeros(R,1);
for i =1:R
    df_r = zeros(R,2); 
    for t=1:2 % 1=backward step; 2=forward step;
        theta_eps = theta0; % "x + epsilon" (small perturbation)
        theta_eps(k) = theta_eps(k) + ((-1)^t)*step;
        [~,PS] = predicted_share(theta_eps,delta0,inp);  % PS is j by r 
        psi = PS(:,i);
        jni = JN(:,i);
        logPSi = log(psi);
        A = logPSi.*jni;
        lklfn = sum(sum(A));
        df_r(i,t) = lklfn;
    end
    fprintf('Current number is %d\n',i);
    dlklfn(i,1) = (df_r(i,2)-df_r(i,1))/(2*step) % diff between second col. and first col.
end