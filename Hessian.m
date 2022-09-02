function [hessian] = Hessian(delta0,theta0,inp)

%%%%% INPUT
dh = inp.dh;
R = inp.R;
nthetap = inp.ntheta_p;
nthetait1 = inp.ntheta_it1;
nthetawd2 = inp.ntheta_wd2;
nthetawd3 = inp.ntheta_wd3;

if inp.ponly
    ntheta =nthetap;
    hessian=zeros(ntheta,ntheta,R);
    
    k=1;
    while k <= ntheta
        jacob = dlkl(delta0,theta0,inp,k); % the sec moment cond. evaluated at optimum theta1 : R by 1
        
        m =1;
        while m <=ntheta
            chg = dh(m);
            theta0(m)=theta0(m)+chg;
            fprintf('Current number is %d\n',k);
            hessian(m,k,:) = (dlkl(delta0,theta0,inp,k)- jacob)./chg
            theta0(m) = theta0(m)-chg;
            m=m+1;
        end
        k = k+1;
    end
elseif inp.inter1
    ntheta=nthetait1;
    hessian=zeros(ntheta,ntheta,R);
    
    k=1;
    while k <= ntheta
        jacob = dlkl(delta0,theta0,inp,k); % the sec moment cond. evaluated at optimum theta1 : R by 1
        
        m =1;
        while m <=ntheta
            chg = dh(m);
            theta0(m)=theta0(m)+chg;
            fprintf('Current number is %d\n',k);
            hessian(m,k,:) = (dlkl(delta0,theta0,inp,k)- jacob)./chg
            theta0(m) = theta0(m)-chg;
            m=m+1;
        end
        k = k+1;
    end
elseif inp.wd2
    ntheta=nthetawd2;
    hessian=zeros(ntheta,ntheta,R);
    
    k=1;
    while k <= ntheta
        jacob = dlkl(delta0,theta0,inp,k); % the sec moment cond. evaluated at optimum theta1 : R by 1
        
        m =1;
        while m <=ntheta
            chg = dh(m);
            theta0(m)=theta0(m)+chg;
            fprintf('Current number is %d\n',k);
            hessian(m,k,:) = (dlkl(delta0,theta0,inp,k)- jacob)./chg
            theta0(m) = theta0(m)-chg;
            m=m+1;
        end
        k = k+1;
    end
    elseif inp.wd3
    ntheta=nthetawd3;
    hessian=zeros(ntheta,ntheta,R);
    
    k=1;
    while k <= ntheta
        jacob = dlkl(delta0,theta0,inp,k); % the sec moment cond. evaluated at optimum theta1 : R by 1
        
        m =1;
        while m <=ntheta
            chg = dh(m);
            theta0(m)=theta0(m)+chg;
            fprintf('Current number is %d\n',k);
            hessian(m,k,:) = (dlkl(delta0,theta0,inp,k)- jacob)./chg
            theta0(m) = theta0(m)-chg;
            m=m+1;
        end
        k = k+1;
    end

end
end