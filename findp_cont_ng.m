function [cpopt,P_opt,dPc_opt,dst] = findp_cont_ng(mc,inp,ownership)
data = inp.d_ng;
dt3 =inp.XS_ng2;

ssp = inp.ssp;
rest = inp.rest_ng;

theta2wd3 = inp.theta2_wd3;
theta1wd3 = inp.theta1_wd3;

XIinter5 = inp.XIinter5;
XIinter5 = XIinter5';


R = inp.R;
J = inp.J;
ninter5 = inp.ninter5;
ninterp4 = inp.ninterp4;



reswd3 = inp.res_ng;

tol = 1.0001e-10;
dst = 1;
cp0 = 1/2*ssp;
cp0([1,4,7]) = 2.5;
if inp.wd3
    theta2= theta2wd3;
    alphac = theta2wd3(2);
    theta1 = theta1wd3;
    res = reswd3;
    while dst > tol
        Exp1 = rest;% drop original price
        Exp2 = [cp0,Exp1];
        Exp3 = [ones(J,1),Exp2]; % total 7 variables
        
        delta_pchg = Exp3*theta2+res; % Jx1
        
        XSp5 = repmat(cp0,1,ninterp4); % Jxninter
        XSinter5= [XSp5,dt3];
        mu = 0;
        for k=1:ninter5 % loop through all product characteristics except price
            mu = mu + theta1(k)*XIinter5(k,:).*repmat(XSinter5(:,k),1,R);
        end
        
        U = repmat(delta_pchg,1,R) + mu;
        eU = exp(U);
        denom = 1+ sum(eU,1);
        P_opt = mean(eU./denom,2);
        PS_opt = eU./denom;
        
        % J x 1 x R
        PP_opt = reshape(PS_opt,J,1,R);
        % J x J x R
        PP_opt = repmat(PP_opt,1,J,1); % Create Nt identical copies of P and stack horizontally
        PPt_opt = permute(PP_opt,[2 1 3]);
        I = repmat(eye(J),1,1,R);
        
        muc = 0;
        for k=1:ninterp4 % loop through all product characteristics except price
            muc = muc + theta1(k)*XIinter5(k,:);
        end
        indiv=muc;
        indiv = reshape(indiv,1,1,R);
        indiv = repmat(indiv,J,1,1);
        indiv = repmat(indiv,1,J,1);
        alpha = repmat(alphac,J,J,R);
        
        coeff = indiv +alpha;
        
        dPc_opt = mean(coeff.*PP_opt.*(I - PPt_opt),3); % NtxNt matrix of derivatives wrt price
        
        mtS_opt = (ownership.*(-dPc_opt))';
        invmtS_opt = inv(mtS_opt);
        OmgP_opt = invmtS_opt*P_opt;
        
        cpopt = cp0 +log(OmgP_opt+mc)-log(cp0)  ;
        cpopt([1,4,7]) = 2.5;
        dst = max(abs(cp0 - cpopt));
        cp0 = cpopt;
    end
    
    
    pctmarkup = (ssp-mc)./ssp;
end