function f = GMM_objective(theta,Z,W,psdiff_fn,inp)
u = psdiff_fn(theta);
R = inp.R;
J = inp.J;

if inp.pinter  %% Use instruments varaince as an optimal weighting matrix
    Zujr = zeros(size(Z,2),1);
    for j = 1:J
        for r = 1:R
            ujr = u(j,r);
            Zj = Z(j,:);
            Zujr = Zujr + Zj'*ujr;
        end
    end
    % for check
    %Zujr0 = zeros(size(Z,2),1);
    %for j = 1:J
        Zujr0 = Z'*sum(u,2);
     %   Zujr0 = Zujr0+Zu;     % the number of instrument by 1
    %end
    Zujr==Zujr0
    
    f = Zujr'*(W\Zujr); % same as Zu'*inv(W)*Zu
 
else                %% Do not use instruments >>>>>>>> Need to be checked
        usqrd=u*u';
        Wfn=sum(sum(usqrd));
        Wfn=Wfn/(J);
        rank(Wfn);
        ujr = sum(sum(u));
        f = ujr'*((1/Wfn)*ujr); 
end

 
