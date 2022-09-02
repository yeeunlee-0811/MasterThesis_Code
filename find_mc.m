function [mc,dPc,pctmarkup]=find_mc(inp,ownership)
data = inp.data;
dt1 = inp.dt1;
dt2 = inp.dt2;

ssp = inp.ssp;
rest = inp.x2rest;


theta2p = inp.theta2_p;
theta1p = inp.theta1_p;
deltap = inp.delta_p;

theta2it1 = inp.theta2_it1;
theta1it1 = inp.theta1_it1;
deltait1 = inp.delta_it1;

theta2wd2 = inp.theta2_wd2;
theta1wd2 = inp.theta1_wd2;
deltawd2 = inp.delta_wd2;

theta2wd3 = inp.theta2_wd3;
theta1wd3 = inp.theta1_wd3;
deltawd3 = inp.delta_wd3;

XIinter1 = inp.XIinter1;
XIinter1 = XIinter1';

XIinter2 = inp.XIinter2;
XIinter2 = XIinter2';

XIinter3 = inp.XIinter3;
XIinter3 = XIinter3';

XIinter4 = inp.XIinter4;
XIinter4 = XIinter4';

XIinter5 = inp.XIinter5;
XIinter5 = XIinter5';

XSinter1 = inp.XSinter1;
XSinter2 = inp.XSinter2;
XSinter3 = inp.XSinter3;
XSinter4 = inp.XSinter4;
XSinter5 = inp.XSinter5;

nup = inp.nup;
nudata = inp.nudata;

R = inp.R;
J = inp.J;
ninter1 = inp.ninter1;
ninter2 = inp.ninter2;
ninter3 = inp.ninter3;
ninter4 = inp.ninter4;
ninter5 = inp.ninter5;
ninterp3 = inp.ninterp3;
ninterp4 = inp.ninterp4;


resp = inp.res_p;
resit1 = inp.res_it1;
reswd2 = inp.res_wd2;


if inp.ponly
    alphac = theta2p(2);
    sigp = theta1p(1);
    sigd = theta1p(2);
    
    delta = deltap; % Jx1
    
    mu = 0;
    mu = mu + sigp*nup.*repmat(ssp,1,R);
    mu = mu + sigd*nudata.*repmat(data,1,R);
    
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = 1+ sum(eU,1);
    P = mean(eU./denom,2);
    PS_opt = eU./denom;
    
    % J x 1 x R
    PP_opt = reshape(PS_opt,J,1,R);
    % J x J x R
    PP_opt = repmat(PP_opt,1,J,1); % Create Nt identical copies of P and stack horizontally
    PPt_opt = permute(PP_opt,[2 1 3]);
    I = repmat(eye(J),1,1,R);
    
    muc = 0;
    muc = muc +sigp.*nup;
    
    indiv=muc;
    indiv = reshape(indiv,1,1,R);
    indiv = repmat(indiv,J,1,1);
    indiv = repmat(indiv,1,J,1);
    alpha = repmat(alphac,J,J,R);
    
    coeff = indiv +alpha;
    
    dPc = mean(coeff.*PP_opt.*(I - PPt_opt),3); % NtxNt matrix of derivatives wrt price
    
    mc = ssp + inv((ownership.*dPc)')*P;
    pctmarkup = (ssp-mc)./ssp;
elseif inp.inter1
    alphac = theta2it1(2);
    theta1 = theta1it1;
    
    delta = deltait1;
    
    mu = 0;
    for k=1:ninter1 % loop through all product characteristics except price
        mu = mu + theta1(k)*XIinter1(k,:).*repmat(XSinter1(:,k),1,R);
    end
    
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = 1+ sum(eU,1);
    P = mean(eU./denom,2);
    PS_opt = eU./denom;
    
    % J x 1 x R
    PP_opt = reshape(PS_opt,J,1,R);
    % J x J x R
    PP_opt = repmat(PP_opt,1,J,1); % Create Nt identical copies of P and stack horizontally
    PPt_opt = permute(PP_opt,[2 1 3]);
    I = repmat(eye(J),1,1,R);
    
    muc = 0;
    for k=1:ninter1 % loop through all product characteristics except price
        muc = muc + theta1(k)*XIinter1(k,:);
    end
    indiv=muc;
    indiv = reshape(indiv,1,1,R);
    indiv = repmat(indiv,J,1,1);
    indiv = repmat(indiv,1,J,1);
    alpha = repmat(alphac,J,J,R);
    
    coeff = indiv +alpha;
    
    dPc = mean(coeff.*PP_opt.*(I - PPt_opt),3); % NtxNt matrix of derivatives wrt price
    mc = ssp + inv((ownership.*dPc)')*P;
    pctmarkup = (ssp-mc)./ssp;
   
elseif inp.wd2
    alphac = theta2wd2(2);
    theta1 = theta1wd2;
    
    delta = deltawd2; % Jx1
    
    mu = 0;
    for k=1:ninter4 % loop through all product characteristics except price
        mu = mu + theta1(k)*XIinter4(k,:).*repmat(XSinter4(:,k),1,R);
    end
    
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = 1+ sum(eU,1);
    P = mean(eU./denom,2);
    PS_opt = eU./denom;
    
    % J x 1 x R
    PP_opt = reshape(PS_opt,J,1,R);
    % J x J x R
    PP_opt = repmat(PP_opt,1,J,1); % Create Nt identical copies of P and stack horizontally
    PPt_opt = permute(PP_opt,[2 1 3]);
    I = repmat(eye(J),1,1,R);
    
    muc = 0;
    for k=1:ninterp4 
        muc = muc + theta1(k)*XIinter4(k,:);
    end
    indiv=muc;
    indiv = reshape(indiv,1,1,R);
    indiv = repmat(indiv,J,1,1);
    indiv = repmat(indiv,1,J,1);
    alpha = repmat(alphac,J,J,R);
    
    coeff = indiv +alpha;
    
    dPc = mean(coeff.*PP_opt.*(I - PPt_opt),3); % NtxNt matrix of derivatives wrt price
    
    mc = ssp + inv((ownership.*dPc)')*P;
    pctmarkup = (ssp-mc)./ssp;
elseif inp.wd3
    alphac = theta2wd3(2);
    theta1 = theta1wd3;
    
    delta = deltawd3; % Jx1
    
    mu = 0;
    for k=1:ninter5 % loop through all product characteristics except price
        mu = mu + theta1(k)*XIinter5(k,:).*repmat(XSinter5(:,k),1,R);
    end
    
    U = repmat(delta,1,R) + mu;
    eU = exp(U);
    denom = 1+ sum(eU,1);
    P = mean(eU./denom,2);
    PS_opt = eU./denom;
    
    % J x 1 x R
    PP_opt = reshape(PS_opt,J,1,R);
    % J x J x R
    PP_opt = repmat(PP_opt,1,J,1); % Create Nt identical copies of P and stack horizontally
    PPt_opt = permute(PP_opt,[2 1 3]);
    I = repmat(eye(J),1,1,R);
    
    muc = 0;
    for k=1:ninterp4
        muc = muc + theta1(k)*XIinter5(k,:);
    end
    indiv=muc;
    indiv = reshape(indiv,1,1,R);
    indiv = repmat(indiv,J,1,1);
    indiv = repmat(indiv,1,J,1);
    alpha = repmat(alphac,J,J,R);
    
    coeff = indiv +alpha;
    
    dPc = mean(coeff.*PP_opt.*(I - PPt_opt),3); % NtxNt matrix of derivatives wrt price
    
    mc = ssp + inv((ownership.*dPc)')*P;
    pctmarkup = (ssp-mc)./ssp;
    
end