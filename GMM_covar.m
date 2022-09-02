function V = GMM_covar(theta,Z,inp,psdiff_fn,heterosked)

J = inp.J;
XS = inp.XS;
K = size(XS,2);
R = inp.R;

[nj,L] = size(Z);
G = nj/J;
u = psdiff_fn(theta);

if heterosked %% slide equation (33) 
    V = zeros(L,L);
    %for k=1:R
    %    ui=u(k);
        for i=1:J
        ix = G*(i-1)+1:G*i;
        temp = Z(ix,:)'*u(ix);
        V = V + temp*temp';
        end
    %end
else
    uu = zeros(G,G);
    %for k = 1:R
    %    ui=u(k);
        for i=1:J
        ix = G*(i-1)+1:G*i;
        temp = u(ix);
        uu = uu + temp*temp';
        end
    %end 
    uu = uu/J; % this is a sample average. Here it is more standard to
    % divide by N instead of N-K, but if we want to check that our code gives
    % the same results as Matlab's inbuilt fitlm when using X as
    % instruments, this requires N-K. Most often, for system estimation,
    % we will not assume homoskedasticity in any case. 
    
    V = zeros(L,L);
    for i=1:J
        ix = G*(i-1)+1:G*i;
        temp = Z(ix,:)';
        V = V + temp*uu*temp';
    end
end
    
