function [coef, stderr, pval,u_iv,PWX] = tsls(y,X,W)
% Original model: y = X*beta + u
% We are concerned that one or more regressors in X
% could be endogenous. 
% W contains the instruments and all regressors in X 
% except the suspected one(s).
% An example: 
% m_t = b1 + b2*r_t + b3*y_t + b4*m_(t-1) + b5*m_(t-2) + u_t
% In this model, we have a concern that r_t could be correlated 
% with u_t and we want to see the results with instruments, 
% r_(t-1) and r_(t-2). So the inputs y, X and W must be
% y = m_t;
% X = [ones(length(y)) r y m_(t-1) m_(t-2)];
% W = [ones(length(y)) r_(t-1) r_(t-2) y m_(t-1) m_(t-2)];
% Please do not forget checking the length of the variables be equal 
% Number of observations
n = length(y);
% Number of the regressors in the original model
k = size(X,2); 
% Degree of freedom
df = n - k; 
% Stage 1 of 2SLS: Regress X on W, 
% and then get fitted values PWX
PWX = W*inv(W'*W)*W'*X;
% Stage 2 of 2SLS: Regress y on fitted values PWX 
beta_2sls = inv(PWX'*PWX)*PWX'*y;
coef      = beta_2sls;
% Structural residuals
u_iv = y - X*beta_2sls;
% Asymptotic covariance
s_hat   = u_iv'*u_iv/df; 
% Estimated covariance matrix
var_hat = s_hat*inv(PWX'*PWX); 
% Standard errors
stderr  = sqrt(diag(var_hat)); 
% t-ratios
t_stat  = beta_2sls./stderr; 
% p-values
pval = betainc(df./(df+(1.*t_stat.^2)),(df./2),(1./2)); 
