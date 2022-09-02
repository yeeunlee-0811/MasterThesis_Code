function [beta,u] = GMM_beta(Y,X,Z,W)
% GMM_beta calculates the gmm estimates (function is called by GMM.m) and
% residuals

P = X'*Z*(W\Z'); % we will use this more than once. A\B = inv(A)*B; either way works.
beta = (P*X)\(P*Y); 
u = Y - X*beta; % residuals