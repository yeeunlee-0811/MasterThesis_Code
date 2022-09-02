function [foc,dP,P]=find_p(price,inp,mc)
% This function generates profit maximizing firms' pricing condition which
% is (Omega')*(p-mc) + S = 0

% Omega : element by element product of I and delta S
% delta S : J by J Predicted share change of product j(row) by price of
% product i(column)
% p : J by 1 price matrix
% S : J by 1 share vector
theta1=inp.theta1;
theta2 = inp.theta;
res2= inp.res2; % Jx1
XSred = inp.XSred; % Explanatory variables without constant (6 variables)
XIinter = inp.XIinter;
XIinter= XIinter';
alpha = theta1(2,:);
J = inp.J;
ninter= inp.ninter;
R = inp.R;
Own=inp.Own;
delta = inp.delta;
S = inp.S;
data =inp.data;
dt = inp.dt;
ssp = inp.ssp;
nup= inp.nup;
nudata = inp.nudata;
XSinter = inp.XSinter;
ninterp=inp.ninterp;

if inp.usedhp
    % XSinter matrix with given price (J by 1)
    Exp1 = data; 
    Exp12 = [price,Exp1];
    Exp13 = [ones(J,1),Exp12]; % total 3 variables
    
    delta_pchg = Exp13*theta1+res2; % Jx1
    
    XSinter1 = repmat(price,1,ninterp); % Jxninter
    XSinterp1 = [XSinter1,dt];
    inp.XSinter = XSinterp1;
    inp.ssp = price;
    [P0] = predicted_share(theta2,delta_pchg,inp);
    

    dhp = stepPrice(inp);
    dP = zeros(J,J);
    k=1;
    while k <= J
        chg = dhp(k);
        price(k) = price(k) + chg;
        % change XSinter matrix
        XSinter2 = repmat(price,1,ninterp);
        XSinterp2 = [XSinter2,dt];
        inp.XSinter = XSinterp2;
        Exp1 = data; 
        Exp1_1 = [price,data];
        Exp1_2 = [ones(J,1),Exp1_1]; % total 7 variables
        
        delta_pchg_2 = Exp1_2*theta1+res2; % Jx1
        inp.ssp = price;
        [Ph] = predicted_share(theta2,delta_pchg_2,inp);
        
        inp.XSinter = XSinter;
        inp.ssp =ssp;
        dP(:,k) = (Ph-P0).*(1/chg);
        
        price(k) =price(k)-chg;
        k = k + 1;
    end
    
    P = P0;
    foc = ((Own.*dP)')*(price-mc) + P;
    
    checkfoc = 1;
else
    thetafin = theta2(1:ninter);
    sigp = theta2(ninter+1);
    sigd = theta2(ninter+2);
    
    Exp1 = data; % drop original price
    Exp122 = [price,Exp1];
    Exp133 = [ones(J,1),Exp122]; % total 7 variables
    
    delta_pchg = Exp133*theta1; % Jx1
    
    XSinter1 = repmat(price,1,ninterp); % Jxninter
    XSinter = [XSinter1,dt];
    mu = 0;
    for k=1:ninter % loop through all product characteristics except price
        mu = mu + thetafin(k)*XIinter(k,:).*repmat(XSinter(:,k),1,R);
    end
    
    mu = mu + sigp*nup.*repmat(ssp,1,R);
    mu = mu + sigd*nudata.*repmat(data,1,R);
    
    U = repmat(delta_pchg,1,R) + mu;
    eU = exp(U);
    denom = 1+sum(eU,1);
    P = mean(eU./denom,2);
    PS = eU./denom;
    
    % J x 1 x R
    PP = reshape(PS,J,1,R);
    % J x J x R
    PP = repmat(PP,1,J,1); % Create Nt identical copies of P and stack horizontally
    PPt = permute(PP,[2 1 3]);
    I = repmat(eye(J),1,1,R);
    
    muc = 0;
    for k=1:ninterp % loop through all product characteristics except price
        muc = muc + thetafin(k)*XIinter(k,:);
    end
    muc = muc +sigp.*nup;
    
    indiv=muc;
    indiv = reshape(indiv,1,1,R);
    indiv = repmat(indiv,J,1,1);
    indiv = repmat(indiv,1,J,1);
    alpha = repmat(alpha,J,J,R);
    coeff = indiv+alpha;
    
    dP = mean(coeff.*PP.*(I - PPt),3); % NtxNt matrix of derivatives wrt price
    
    foc = ((Own.*dP)')*(price-mc) + P; % + penalty;
end