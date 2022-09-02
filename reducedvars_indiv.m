clear
%% Importing data
% Importing individual data
panelm = importdata('exp_est1.xls');
XId = panelm;
indiv_varnames = XId(1,:);
XId = array2table(XId, 'VariableNames',indiv_varnames);
XId(:,1)=[]; % drop person ID
XId(1,:)=[]; % drop varnames

[R,nikrc]=size(XId);

%% Importing Service data
% In the reduced version, all products are classifed into one of nine
% products, mostly classified by their data and by their providers
XSd = importdata('spconv.xls');
XS = XSd.data;
XSvarnames = XSd.textdata(1,:) ;
XSvarnames(:,1)=[];
XSd = array2table(XS, 'VariableNames', XSvarnames);

%% <<<<<<<<<<<<<<<<< Reducing products >>>>>>>>>>>>>>>>>>>
%%% First, adjust products
schoice = XId(:,{'schoice15'});
schoice=table2array(schoice);
schoice=str2double(schoice);
uniqck1 =unique(schoice);
% Adjusting schoice after removing the products whose shares are zero
for i = 1:2455
    for d1 = 31:34
        if schoice(i,1) ==d1
           schoice(i,1) = d1-2;
        end
    end
    for d2=36:41
        if schoice(i,1)==d2
            schoice(i,1) = d2-3;
        end
    end
    for d3= 43:58
        if schoice(i,1)==d3
            schoice(i,1)=d3-4;
        end
    end
    if schoice(i,1) ==60
        schoice(i,1)=55;
    end
    for d4 = 62:63
        if schoice(i,1)==d4
            schoice(i,1) = d4-6;
        end
    end
end
uniqck2 = unique(schoice);

% Change to 9 products
redschoice=zeros(R,1);
for i =1:R
    for d1=[1:3,7:9,14]
        if schoice(i,1) ==d1
            redschoice(i,1) =1;
        end
    end
     for d2=[4:6,10:13]
        if schoice(i,1) ==d2
            redschoice(i,1) =2;
        end
     end
    for d3=15:17
        if schoice(i,1) ==d3
            redschoice(i,1) =3;
        end
    end
    for d4 = [23:25,33:35]
        if schoice(i,1)==d4
            redschoice(i,1)=4;
        end
    end
    for d5 = [18:21,26:28,29,36]
        if schoice(i,1)==d5
            redschoice(i,1) =5;
        end
    end
    for d6 = [22,30:32,37:38]
        if schoice(i,1)==d6
            redschoice(i,1)=6;
        end
    end
    for d7 = 43:51
        if schoice(i,1)==d7
            redschoice(i,1)=7;
        end
    end
    for d8 =[39:42,52:54]
        if schoice(i,1)==d8
            redschoice(i,1)=8;
        end
    end
    for d9 =55:57
        if schoice(i,1)==d9
            redschoice(i,1)=9;
        end
    end
end
uniqck3 = unique(redschoice);
isempty(redschoice(:))

%% Setting products data
XS = table2array(XSd);
row = find(XS(:,19)==0);
XS(row,:) =[];
S=XS(:,19);
s_red=XS(:,23);
share_red=zeros(9,1);
for i =1:9
    same_sred = s_red ==i;
    diff_sred = s_red ~=i;
    share_red(i,1)=sum(S(same_sred));
end
all(share_red>0)

S = share_red;
logS =log(share_red);

% turning to reduced products data
XStored = XS;
% case1 : data, intrafree, interfree, price firm dummy (skt, kt): total 6
% --> In this case, price's pvalue = 0.6~
XStored(:,[2,5:17,19,22:end])=[];

% case2 : data, intrafree, interfree, price firm dummy (skt, kt)
%XStored(:,[2,8,10:17,19,22:end])=[];

% service dataset under 9 products case 
XSred=zeros(9,size(XStored,2));
for i=1:9
     same_sred = s_red==i ;
     diff_sred = s_red~=i ;
    for l=1:size(XStored,2)
        zk = XStored(:,l);
        XSred(i,l) = mean(zk(same_sred)) ;  
    end
end
% adjust DATA unit to solve inf problem
XSred(:,1)=XSred(:,1);  %% GB to Terrabyte
[J,nskrc]=size(XSred);

%% Interaction terms
% Income dummy
income = XId(:,{'inc15'});
income = table2array(income);
income = str2double(income);

dincome =zeros(size(XId,1),3); 
for k=1:3
    di = income == k+1;
    dincome(:,k)=di;
end

payer = XId(:,{'payer1'});
payer = table2array(payer);
payer = str2double(payer);
incpay = dincome.*payer; 

ind = payer==0;
npayer = ind;
% age dummy
age = XId(:,{'age15'});
age = table2array(age);
age = str2double(age);

dage = zeros(size(XId,1),4);
for k=1:4
    ida = age == k;
    dage(:,k)=ida;
end
agenpay = dage.*npayer;

XSred(:,1) = XSred(:,1)./10;
data = XSred(:,1);
intra= XSred(:,2);
inter= XSred(:,3);
inp.data = data;
inp.intra = intra;
inp.inter = inter;

Xit = XId(:,{'hhldsiz15', 'gender', 'mar2','djob'});
Xit = table2array(Xit);
Xit = str2double(Xit);

Xindiv = Xit;
nindiv = size(Xindiv,2);
ninterp = size([incpay,agenpay],2);
XIinter = [incpay,agenpay,Xindiv,Xindiv,Xindiv]; % R by 3+4

XSinter1= repmat(XSred(:,4),1,ninterp); %  J(9) by 3+4
dt = repmat(data,1,nindiv);
dntra =repmat(intra,1,nindiv);
dnter = repmat(inter,1,nindiv);

XSinter=[XSinter1,dt,dntra,dnter];


ninter =size(XIinter,2);

inp.nindiv=nindiv;
inp.ninterp=ninterp;
inp.ninter=ninter;

JN = zeros(J,R);
 for j = 1:J
     for r=1:R
         if redschoice(r,1)==j
         JN(j,r) = 1;
         end
     end
 end
%% Input 
inp.JN =JN;
inp.R =R;
inp.J =J;
inp.XIinter= XIinter;
inp.XSinter= XSinter;
inp.XS =XSred;
inp.nskrc=nskrc;
inp.ninter=ninter;
inp.nikrc=nikrc;
inp.logS = logS;
inp.S = S;
%% Ownership matrix 

% Ownership matrix : J BY J
Ipre = [1;1;1;2;2;2;3;3;3];
Own = zeros(J,J);
for i =1:J
    Own(:,i) = Ipre ==Ipre(i);
end
inp.Own = Own;

%% <<<<<<<<<<<<<<<<<<<<<<<<< Likelilhood >>>>>>>>>>>>>>>>>>>>>>>>>>>
%% 
obj = @(theta) lklfn2_indiv(inp,theta);
options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

theta0_indiv = 0.5*ones(ninter,1);
thetafin_indiv = fminunc(obj, theta0_indiv, options);


[~,~,~,deltafin_indiv] = lklfn2(inp,thetafin_indiv);
[Pfin_indiv] = predicted_share(thetafin_indiv,deltafin_indiv,inp);

ntheta1_indiv = size(thetafin_indiv,1); % ninter + nikrc
parmnum1_indiv = ntheta1_indiv + J;
% Explanatory variables
%% 
spbd = XS(:,10)/10000;
spbd_red=zeros(9,1);
for i=1:9
     same_sred = s_red==i ;
     spbd_red(i,1) = mean(spbd(same_sred)) ;  
end
%changing the order of variables
ssp =XSred(:,4);
XSred0=XSred;
XSred0(:,4)=[];
XSred=[ssp, XSred0];
% Explanatory variables
%%
X=[ones(J,1),XSred];
INS = [ones(J,1),spbd_red,XSred0]; %----------------------------------------> Need to check IV !!!

[theta1_indiv,stderr_indiv,pval_indiv,u_iv_indiv]=tsls(deltafin_indiv, X,INS);

%% Inputr Structure

inp.theta1= theta1_indiv;
inp.delta = deltafin_indiv;
inp.theta = thetafin_indiv;

inp.X = X;
inp.XSred = XSred;
inp.res2 = u_iv_indiv;
inp.Pshare = Pfin_indiv;
inp.eps = 1.0e-6;
inp.ntheta = ntheta1_indiv;
inp.parmnum1 = parmnum1_indiv;
inp.Dlikelmat0 = zeros(R,ntheta1_indiv);
inp.ssp = ssp;

[dh] = stepP(inp);
inp.dh = dh;

[dhd] = stepD(inp);
inp.dhd = dhd;

%% Calculating se
[V1_indiv, Dlikelmat_indiv]=varMOM(inp);
%% 
Gam_indiv = Gamma(inp, Dlikelmat_indiv);
%% 

G2_indiv = Gam_indiv(J+1:J+ntheta1_indiv, J+1:J+ntheta1_indiv);
V4_indiv = V1_indiv(J+1:J+ntheta1_indiv,J+1:J+ntheta1_indiv);
VW_indiv = inv(G2_indiv'*G2_indiv)*G2_indiv'*V4_indiv*G2_indiv*inv(G2_indiv'*G2_indiv)./R;
sew_indiv = diag(sqrt(VW_indiv));


%% price according to foc rule with zero mc
inp.dt=dt;
inp.usedhp = true;
mc =0;
findp = @(price) find_p(price,inp,mc);

 options = optimoptions('fsolve','Display','iter',...
     'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
     'StepTolerance',1e-12,'FunctionTolerance',1e-12,...
     'MaxFunctionEvaluations', 1e+12,'MaxIterations', 4e+4 );

price0 = ssp.*(1/3);

optprice = fsolve(findp,price0,options);

[~,~,optP]=find_p(optprice,inp,0);
diffP = optP - S;

%% implied markups / marginal cost
%find mc with assump. that firm set price based on foc rule

inp.usedhp = true;
findmc = @(mc) find_p(ssp,inp,mc);
mc0 = ssp.*(1/2);

 options = optimoptions('fsolve','Display','iter',...
     'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
     'StepTolerance',1e-12,'FunctionTolerance',1e-12,...
     'MaxFunctionEvaluations', 1e+12,'MaxIterations', 4e+4 );
 
mccal = fsolve(findmc,mc0,options);

%% Elasticities 
inp.usedhp=true;
mc=0; % It won't be used to calculate ela though. 
[~,Ela] = find_p(ssp,inp,mc);
[~,Ela2] = pricederiv(ssp,inp);

orgshare = repmat(S,1,J); % replicate Jx1
orgprice = repmat(ssp',J,1); % replicate 1xJ 
orgshare = 1./orgshare;
orgPSshare = orgprice.*orgshare;

elasticity = Ela.*orgPSshare;
