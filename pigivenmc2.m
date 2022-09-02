function [pi2] = pigivenmc2(mc2,inp,ownership)

ssp = inp.ssp;
res2 = inp.res2;
data = inp.data;
thetafin_indiv = inp.theta;
theta1_indiv = inp.theta1;
ninterp = inp.ninterp;
XIinter = inp.XIinter;
XIinter = XIinter';
nudata = inp.nudata;
nup = inp.nup;
R = inp.R;
J = inp.J;
ninter = inp.ninter;
dt = inp.dt;
Own = inp.Own;

tol = 1.0001e-10;
dst = 1;

cprice0_opt = 1/2*ssp; 
deltac = zeros(J,1);
alphac = theta1_indiv(2,:);
theta2 = thetafin_indiv(1:ninter);
sigp = thetafin_indiv(ninter+1);
sigd = thetafin_indiv(ninter+2);


while dst > tol
    Exp1 = data; % drop original price
    Exp122 = [cprice0_opt,Exp1];
    Exp133 = [ones(J,1),Exp122]; % total 7 variables
    
    delta_pchg = Exp133*theta1_indiv+res2; % Jx1
    
    XSinter1 = repmat(cprice0_opt,1,ninterp); % Jxninter
    XSinter = [XSinter1,dt];
    mu = 0;
    for k=1:ninter % loop through all product characteristics except price
        mu = mu + theta2(k)*XIinter(k,:).*repmat(XSinter(:,k),1,R);
    end
    
    mu = mu + sigp*nup.*repmat(ssp,1,R);
    mu = mu + sigd*nudata.*repmat(data,1,R);
    
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
    for k=1:ninterp % loop through all product characteristics except price
        muc = muc + theta2(k)*XIinter(k,:);
    end
    muc = muc +sigp.*nup;
    
    indiv=muc;
    indiv = reshape(indiv,1,1,R);
    indiv = repmat(indiv,J,1,1);
    indiv = repmat(indiv,1,J,1);
    alpha = repmat(alphac,J,J,R);
    
    coeff = indiv +alpha;
    
    dPc_opt = mean(coeff.*PP_opt.*(I - PPt_opt),3); % NtxNt matrix of derivatives wrt price
    
    Ownmono = ones(9,9);
    
    mtS_opt = (ownership.*(-dPc_opt))';
    invmtS_opt = inv(mtS_opt);
    OmgP_opt = invmtS_opt*P_opt;
    
    % <<<<<<<<<<<<<<<<< SET MC >>>>>>>>>>>>>>>>>>>>>>
    mc1 = inp.mc1;
    % mc2 is given 
    optp3 = @(mc3) pigivenmc3(mc3,mc2,inp,ownership);
    
    options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

    Aforineq3 = -eye(3);
    bforineq3 = zeros(3,1);
    Aforeqconst3 = ones(1,3);
    bforeqconst3 = inp.eqconst3;
    mc30= 4*ones(3,1);
    optmccal3= fmincon(optp3,mc30,Aforineq3,bforineq3,Aforeqconst3,bforeqconst3,[],[],[]);
    
    % <<<<<<<<<<<<<<< MC after firm 3's optimization >>>>>>>>>>>>>>>>
    mc = [mc1;mc2;optmccal3];
    
    cprice_opt = cprice0_opt +log(OmgP_opt+mc)-log(cprice0_opt)  ;
    dst = max(abs(cprice0_opt - cprice_opt));
    cprice0_opt = cprice_opt;
end

    piperj = R.*P_opt.*(cprice_opt - mc);
    pi1 = -sum(piperj(1:3));
    pi2 = -sum(piperj(4:6));
    pi3 = -sum(piperj(7:9));
end