function stdErr = GMM_se(theta,Z,inp,psdiff_fn,linear,heterosked)

V = GMM_covar(theta,Z,inp,psdiff_fn,heterosked);
if linear % if model has the form Y = X*theta + u
    X = inp.XS;
    ZX = Z'*X;
    % varTheta = inv(ZX'*inv(V)*ZX);
    varTheta = inv(ZX'*(inv(Z)*X)); % to avoid Matlab's warning about inverses...
else % Use finite-difference derivatives 
    % Calculate LxK matrix of derivatives of Z'u with respect to each parameter:
    step = 1e-5; % step-length for forward-backward finite differences
    K = length(theta); 
    L = size(Z,2);
    J = inp.J;
    DZu = zeros(L,K); % This 
    for k=1:K
        disp(['Calc. deriv. wrt. parameter ' num2str(k) ' of ' num2str(K)])
        DZu_k = zeros(L,2);
        for t=1:2 % 1=backward step; 2=forward step; 
            theta_eps = theta; % "theta + epsilon" (small perturbation)
            theta_eps(k) = theta_eps(k) + ((-1)^t)*step; 
            DZu_k(:,t) = Z'*psdiff_fn(theta_eps); 
        end
        DZu(:,k) = diff(DZu_k,1,2)/(2*step); % diff between second col. and first col. 
    end
    % A consistent estimator. Same expression as in Wooldridge: Econometric 
    % Analysis of Cross Section and Panel Data, 2nd ed., eq. (14.18), but
    % with V evaluated at the final estimates rather than preliminary
    % estimates (either way gives a consistent estimator), and letting N's 
    % cancel:
    % varTheta = inv(DZu'*inv(V)*DZu);
    all(eig(V)>0)
    all(eig(inv(V))>0)
    ggg =DZu'*(inv(V)*DZu);
    rank(DZu)
    interim = DZu'*(V\DZu);
    issymmetric(interim)
    interimeig =eig(interim)
    varTheta = inv(DZu'*(inv(V)*DZu));
end
stdErr = sqrt(diag(varTheta));
