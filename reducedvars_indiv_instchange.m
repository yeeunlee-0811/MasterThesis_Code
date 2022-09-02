clear
%% Importing data
% Importing individual data
panelm = importdata('exp_est1.xls');
panelmt = importdata('exp_est1.xls');
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

%%
% Change to 9 products
redschoice=zeros(R,1);
for i =1:R
    for d1=[1:3,7:9]
        if schoice(i,1) ==d1
            redschoice(i,1) =1;
        end
    end
     for d2=[4:6,10:11]
        if schoice(i,1) ==d2
            redschoice(i,1) =2;
        end
     end
    for d3=[12:13,15:17]
        if schoice(i,1) ==d3
            redschoice(i,1) =3;
        end
    end
    for d4 = [23:25,33:35]
        if schoice(i,1)==d4
            redschoice(i,1)=4;
        end
    end
    for d5 = [20:21,26:27,29]
        if schoice(i,1)==d5
            redschoice(i,1) =5;
        end
    end
    for d6 = [18:19,22,28,30:32,36,37:38]
        if schoice(i,1)==d6
            redschoice(i,1)=6;
        end
    end
    for d7 = 43:51
        if schoice(i,1)==d7
            redschoice(i,1)=7;
        end
    end
    for d8 =[39:40,52:53]
        if schoice(i,1)==d8
            redschoice(i,1)=8;
        end
    end
    for d9 =[41:42,54,55:57]
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
XS(row,:) = [];
S=XS(:,19);
s_red=XS(:,23);
s_red(14,:)
s_red(14,:) = 0;
share_red=zeros(9,1);
%% 
for i =1:9
    same_sred = s_red ==i;
    diff_sred = s_red ~=i;
    share_red(i,1)=sum(S(same_sred));
end
all(share_red>0)

S = share_red; % share of 9 products except outside good
logS =log(share_red);
%% 
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
inp.income = income;

payer = XId(:,{'payer1'});
payer = table2array(payer);
payer = str2double(payer);
incpay = dincome.*payer; 
incxpay = income.*payer;
inp.incxpay = incxpay;

ind = payer==0;
npayer = ind;

% age dummy
age = XId(:,{'age15'});
age = table2array(age);
age = str2double(age);

agexnpay = age.*npayer;
inp.agexnpay = agexnpay;

dage = zeros(size(XId,1),4);
for k=1:4
    ida = age == k;
    dage(:,k)=ida;
end
agenpay = dage.*npayer;

djob = XId(:,{'djob'});
djob = table2array(djob);
djob = str2double(djob);

edu = XId(:,{'school15'});
edu = table2array(edu);
edu = str2double(edu);

inp.age = age;
inp.djob = djob;
inp.edu = edu;

XSred(:,1) = XSred(:,1)./10;
ssp = XSred(:,4);
data  = XSred(:,1); 
intra= XSred(:,2);
inter= XSred(:,3);
inp.data = data;
inp.intra = intra;
inp.inter = inter;
inp.ssp = ssp;
dt1 = repmat(data,1,3);
inp.dt1=dt1;
dt2 = repmat(data,1,3);
inp.dt2=dt2;
dt3 = repmat(data,1,2);
inp.dt3 = dt3;

dda = dage(:,2:end);
XIinter_wodata1 = incpay; 
XIinter_wodata2 = [incpay,agenpay];
XIinter_wdata1=[incpay,dda];
XIinter_wdata2 = [incpay,age,edu,djob];
XIinter_wdata3 = [incpay,age,djob];
XIinter_wdata4 = [incpay,age,edu];

ninter1 = size(incpay,2);
ninter2 = size([incpay,agenpay],2);
ninterp3 = ninter1;
ninterp4 = ninter1;
ninter3 = size(XIinter_wdata1,2);
ninter4 = size(XIinter_wdata2,2);
ninter5 = size(XIinter_wdata3,2);

XSinter1 = repmat(XSred(:,4),1,size(incpay,2));
XSinter2= repmat(XSred(:,4),1,ninter2); %  J(9) by 3+4
XSinter3 = [repmat(ssp,1,ninterp3),dt1];
XSinter4 = [repmat(ssp,1,ninterp3),dt2];
XSinter5 = [repmat(ssp,1,ninterp3),dt3];
XSinter6 = [repmat(ssp,1,ninterp3),dt3];

XSinter_wodata1= XSinter1;
XSinter_wodata2 = XSinter2;

inp.ninter1 = ninter1;
inp.ninter2 = ninter2;
inp.ninter3 = ninter3;
inp.ninter4 = ninter4;
inp.ninter5= ninter5;
inp.ninterp3 = ninterp3;
inp.ninterp4 = ninterp4;

inp.XSinter1 = XSinter1;
inp.XSinter2 = XSinter2;
inp.XSinter3 = XSinter3;
inp.XSinter4 = XSinter4;
inp.XSinter5= XSinter5;
inp.XSinter6 = XSinter6;

inp.XIinter1 = XIinter_wodata1;
inp.XIinter2 = XIinter_wodata2;
inp.XIinter3 = XIinter_wdata1;
inp.XIinter4 = XIinter_wdata2;
inp.XIinter5 = XIinter_wdata3;
inp.XIinter6 = XIinter_wdata4;

JN = zeros(J,R);
 for j = 1:J
     for r=1:R
         if redschoice(r,1)==j
         JN(j,r) = 1;
         end
     end
 end
 
%%
XSred0=XSred;
XSred0(:,4)=[];
XSred=[ssp, XSred0];

rng(10);
rcp = randn(1,R);
nup = exp(rcp);
rng(11);
rcdata = randn(1,R);
nudata = rcdata;
%% Input 

inp.nup =nup;
inp.nudata= nudata;

inp.XSred0 = XSred0;
inp.ssp = ssp;
inp.JN =JN;
inp.R =R;
inp.J =J;
inp.XS =XSred;
inp.nskrc=nskrc;
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
%% 1. random price & random data only

inp.ponly =1;
inp.inter1=0;
inp.wd2=0;

obj2 = @(theta) lklfn2_ponly(inp,theta);
options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

theta0_ponly = 0.5*ones(2,1);
thetafin_ponly = fminunc(obj2, theta0_ponly, options);

[~,~,~,deltafin_ponly] = lklfn2_ponly(inp,thetafin_ponly);
[Pfin_ponly] = predicted_share(thetafin_ponly,deltafin_ponly,inp);

ntheta1_ponly = size(thetafin_ponly,1); % ninter + nikrc
parmnum1_ponly = ntheta1_ponly + J;

%% 4. INTER1

inp.ponly =0;
inp.inter1=1;
inp.wd2=0;

obj4 = @(theta) lklfn2_inter(inp,theta);
options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

theta0_it1 = 0.5*ones(ninter1,1);
thetafin_it1 = fminunc(obj4, theta0_it1, options);

[~,~,~,deltafin_it1] = lklfn2_inter(inp,thetafin_it1);
[Pfin_it1] = predicted_share(thetafin_it1,deltafin_it1,inp);

ntheta1_it1 = size(thetafin_it1,1); % ninter + nikrc
parmnum1_it1 = ntheta1_it1 + J;

%% 7. Wdata1 : inc*pay + data*(age, edu, djob)

inp.ponly =0;
inp.inter1=0;
inp.wd2 = 1;

obj7 = @(theta) lklfn2_wd(inp,theta);
options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

theta0_wd2 = 0.5*ones(ninter4,1);
thetafin_wd2 = fminunc(obj7, theta0_wd2, options);

[~,~,~,deltafin_wd2] = lklfn2_wd(inp,thetafin_wd2);
[Pfin_wd2] = predicted_share(thetafin_wd2,deltafin_wd2,inp);

ntheta1_wd2 = size(thetafin_wd2,1); % ninter + nikrc
parmnum1_wd2 = ntheta1_wd2 + J;
%% 7. Wdata1 : inc*pay + data*(age, djob)

inp.ponly =0;
inp.inter1=0;
inp.wd2 = 0;
inp.wd3 = 1;
inp.wd4 =0;

obj7 = @(theta) lklfn2_wd(inp,theta);
options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-12,'StepTolerance',1e-12);

theta0_wd3 = 0.5*ones(ninter5,1);
thetafin_wd3 = fminunc(obj7, theta0_wd3, options);

[~,~,~,deltafin_wd3] = lklfn2_wd(inp,thetafin_wd3);
[Pfin_wd3] = predicted_share(thetafin_wd3,deltafin_wd3,inp);

ntheta1_wd3 = size(thetafin_wd3,1); % ninter + nikrc
parmnum1_wd3 = ntheta1_wd3 + J;

%% 
% spbd = XS(:,10)/10000;
% spbd_red=zeros(9,1);
% for i=1:9
%      same_sred = s_red==i ;
%      spbd_red(i,1) = mean(spbd(same_sred)) ;  
% end
%changing the order of variables

spay = XId(:,{'mtotpay14'});
spay_int = XId(:,{'mtotpay_intv14'});
spay = table2array(spay);
spay = str2double(spay);
spay_int=table2array(spay_int);
spay_int=str2double(spay_int);
spaynan = isnan(spay(:));
rowspaynan = find(spaynan ==1);
spay(rowspaynan) = [];
spay_int(rowspaynan) =[];
spayredschoice = redschoice;
spayredschoice(rowspaynan)= [];
lagspay = zeros(9,1);
for k = 1:9
    same_s = spayredschoice == k;
    lagspay(k,1) = mean(spay(same_s));
    histogram(spay_int(same_s))
end
lagspay = lagspay./10;
%%
fmid = XSred(:,end-1:end);
X=[ones(J,1),ssp,data, fmid];
INS = [ones(J,1),lagspay,data, fmid]; %----------------------------------------> Need to check IV !!!
x2rest = [data,fmid];
inp.x2rest=x2rest;
%% Change delta and theta

[theta1_ponly,stderr_ponly,pval_ponly,u_iv_ponly]=tsls(deltafin_ponly, X,INS);
[theta1_it1,stderr_it1,pval_it1,u_iv_it1]=tsls(deltafin_it1, X,INS);
[theta1_wd2,stderr_wd2,pval_wd2,u_iv_wd2]=tsls(deltafin_wd2, X,INS);
[theta1_wd3,stderr_wd3,pval_wd3,u_iv_wd3]=tsls(deltafin_wd3, X,INS);

%%

inp.theta2_p = theta1_ponly;
inp.delta_p = deltafin_ponly;
inp.theta1_p = thetafin_ponly;

inp.theta2_it1= theta1_it1;
inp.delta_it1 = deltafin_it1;
inp.theta1_it1 = thetafin_it1;

inp.theta2_wd2= theta1_wd2;
inp.delta_wd2 = deltafin_wd2;
inp.theta1_wd2 = thetafin_wd2;

inp.theta2_wd3= theta1_wd3;
inp.delta_wd3 = deltafin_wd3;
inp.theta1_wd3 = thetafin_wd3;


inp.X = X;
inp.XSred = XSred;

inp.res_p = u_iv_ponly;
inp.res_it1 = u_iv_it1;
inp.res_wd2 = u_iv_wd2;
inp.res_wd3 = u_iv_wd3;

inp.eps = 1.0e-6;

inp.ntheta_p = ntheta1_ponly;
inp.ntheta_it1 = ntheta1_it1;
inp.ntheta_wd2 = ntheta1_wd2;
inp.ntheta_wd3 = ntheta1_wd3;

%inp.Dlikelmat0 = zeros(R,ntheta1);
inp.ssp = ssp;

 %% Variance with Hessian matrix
inp.ponly =1;
inp.inter1=0;
inp.wd2 = 0;

[dh] = stepP(inp);
inp.dh = dh;

[dhd] = stepD(inp);
inp.dhd = dhd;

delta0 = inp.delta_p;
theta0 = inp.theta1_p;
Hp = Hessian(delta0,theta0,inp);
mHp = mean(Hp,3);
Infomp = -mHp;
var_p = inv(Infomp);

se_p = sqrt(diag(var_p))./sqrt(R);
%%
inp.ponly =0;
inp.inter1=1;
inp.wd2 = 0;

[dh] = stepP(inp);
inp.dh = dh;

[dhd] = stepD(inp);
inp.dhd = dhd;

delta0 = inp.delta_it1;
theta0 = inp.theta1_it1;
Hit1 = Hessian(delta0,theta0,inp);
mHit1 = mean(Hit1,3);
Infomit1 = -mHit1;
var_it1 = inv(Infomit1);

se_it1 = sqrt(diag(var_it1))./sqrt(R);
%%
inp.ponly =0;
inp.inter1=0;
inp.wd2 = 1;

[dh] = stepP(inp);
inp.dh = dh;

[dhd] = stepD(inp);
inp.dhd = dhd;

delta0 = inp.delta_wd2;
theta0 = inp.theta1_wd2;
Hwd2 = Hessian(delta0,theta0,inp);
mHwd2 = mean(Hwd2,3);
Infomwd2 = -mHwd2;
var_wd2 = inv(Infomwd2);

se_wd2 = sqrt(diag(var_wd2))./sqrt(R);

%%
inp.ponly =0;
inp.inter1=0;
inp.wd2 = 0;
inp.wd3 = 1;
inp.wd4 =0;

[dh] = stepP(inp);
inp.dh = dh;

[dhd] = stepD(inp);
inp.dhd = dhd;

delta0 = inp.delta_wd3;
theta0 = inp.theta1_wd3;
Hwd3 = Hessian(delta0,theta0,inp);
mHwd3 = mean(Hwd3,3);
Infomwd3 = -mHwd3;
var_wd3 = inv(Infomwd3);

se_wd3 = sqrt(diag(var_wd3))./sqrt(R);


% % Calculating se
% [V1_indiv, Dlikelmat_indiv]=varMOM(inp);
% % 
% Gam_indiv = Gamma(inp, Dlikelmat_indiv);
% 
% G2_indiv = Gam_indiv(J+1:J+ntheta1_indiv, J+1:J+ntheta1_indiv);
% V4_indiv = V1_indiv(J+1:J+ntheta1_indiv,J+1:J+ntheta1_indiv);
% VW_indiv = inv(G2_indiv'*G2_indiv)*G2_indiv'*V4_indiv*G2_indiv*inv(G2_indiv'*G2_indiv)./R;
% sew_indiv = diag(sqrt(VW_indiv));


%% Finding optimal price with contraction

inp.ponly =0;
inp.inter1=0;
inp.wd2 = 0;
inp.wd3 = 1;

ownership = Own; %Own; % ones(9,9)

% Choose one of MC : MC or MC_ASS
mc = zeros(9,1);
% mc is assumed from unit marginal cost of data
mc_ass = [1.0812;2.1120;3.4580;1.03;2.3100;3.4790;1.1125;2.3925;3.6167];


[cpoptwd3_ass,Poptwd3_ass,dPcwd3_ass] = findp_cont(mc_ass,inp,ownership);

pctmarkup_cf = (cpoptwd3_ass-mc_ass)./cpoptwd3_ass;


%% implied markups / marginal cost
%find mc with assump. that firm set price based on foc rule
inp.ponly =0;
inp.inter1=0;
inp.wd2 = 0;
inp.wd3 =1;
ownership = Own ; % Own %ones(9,9);

[mc_wd3,dPc_wd3,pctmarkup_wd3]=find_mc(inp,ownership);

ppmc1= mc_wd3(1:3).*(S(1:3)./sum(S(1:3)));
ppmc2= mc_wd3(4:6).*(S(4:6)./sum(S(4:6)));
ppmc3= mc_wd3(7:9).*(S(7:9)./sum(S(7:9)));
ppmcf1 = sum(ppmc1);
ppmcf2 = sum(ppmc2);
ppmcf3 = sum(ppmc3);

npp = R.*S;
nppF = zeros(3,1);
nppF(1) = sum(npp(1:3));
nppF(2) = sum(npp(4:6));
nppF(3) = sum(npp(7:9));
%%
totmc = zeros(3,1);
totmc(1) = 2500%ppmcf1*nppF(1);
totmc(2) = 2000%ppmcf2*nppF(2);
totmc(3) = 1700%ppmcf3*nppF(3);
inp.totmc = totmc;

%% Calculating optimal mc
inp.ponly =0;
inp.inter1=0;
inp.wd2 = 0;
inp.wd3 = 1;

ownership = Own;
discrate = [0.6,0.9];
fnoptmc=@(unitmc) find_optmc(unitmc,discrate,ownership,inp);

initmc = ones(3,1); % per firm unit mc 
options = optimoptions('fsolve','Display','iter',...
    'SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-12,...
    'StepTolerance',1e-12,'FunctionTolerance',1e-12,...
    'MaxIterations',2000000000000);

optmccal=fsolve(fnoptmc,initmc,options);

[~,optmcspf0,optP_opt0,optpi,cpp0]=find_optmc(optmccal,discrate,ownership,inp);
%%
optmcspf0 = optmcspf0';

%%
[cpoptwd3_ctf0,Poptwd3_cft0,dPcwd3_ctf0] = findp_cont(optmcspf0,inp,ownership);

%% Elasticities  ----> ADJUST dPc

poverS = repmat(ssp',J,1)./repmat(S,1,J);
elasticity = dPc_wd2.*poverS;

%% New good 1 AND RUN WITH WD3
spec_new1 = [2.5, 0.1]; %0.144166666269302];
spec_new2 = [2.5, 0.1];%0.158333333333333];
spec_new3 = [2.5, 0.1];%0.1483333296246];
fmid_ng= [1,0;0,1;0,0];
spec_new =[spec_new1;spec_new2;spec_new3];
XSng_b = [spec_new,fmid_ng];
XSng0 = [ssp,data,fmid];
XSng0(1,:) = XSng_b(1,:);
XSng0(4,:) = XSng_b(2,:);
XSng0(7,:) = XSng_b(3,:);
rest_ng = XSng0(:,2:end);
inp.rest_ng = rest_ng;

p_ng = XSng0(:,1);
d_ng = XSng0(:,2);
inp.p_ng = p_ng;
inp.d_ng = d_ng;
XS_ng1= repmat(p_ng,1,ninterp4); %  J(9) by 3+4
XS_ng2 = repmat(d_ng,1,ninter5-ninterp4);
inp.XS_ng2= XS_ng2;
XS_ng = [XS_ng1,XS_ng2];

theta1wd3= inp.theta1_wd3;
theta2wd3 = inp.theta2_wd3;

XIap = inp.XIinter5; % R by 3+4
XIap = XIap';

% mu_ng = 0;
% for k=1:ninter5 % loop through all product characteristics except price
%     mu_ng = mu_ng + theta1wd3(k)*XIap(k,:).*repmat(XS_ng(:,k),1,R);
% end
% 
% X_ng = [ones(size(XSng0,1),1),XSng0]; %
res0_ng = inp.res_wd3;
res_ng = res0_ng;
inp.res_ng = res_ng;
% delta_ng = X_ng*theta2wd3+res_ng;
% 
% U_ng = repmat(delta_ng,1,R) + mu_ng;
% eU_ng = exp(U_ng);
% denom_ng = 1+sum(eU_ng,1);
% P_ng = mean(eU_ng./denom_ng,2);
% PS_ng = eU_ng./denom_ng;


%% Calculating prices
inp.ponly =0;
inp.inter1=0;
inp.wd2 = 0;
inp.wd3 = 1;
%%
ownership = Own;

ownm1 = zeros(9);
ownm1(1:3,1:3) = 1;
ownm1(4:9,4:9) =1;

ownm2 = ones(9);

mc_im_ng =[0.931906998;3.22739274;4.935778355;1.429208289;4.266644302;5.253710523;1.299566144;4.004606274;6.119289268];% mc_wd3;%[0.931906998;3.22739274;4.935778355;1.429208289;4.266644302;5.253710523;1.299566144;4.004606274;6.119289268];


mc_px_ng = [0.75;2.112;3.458;0.65;2.31;3.479;0.75;2.3925;3.616666667];%mc_ass;%[0.75;2.112;3.458;0.65;2.31;3.479;0.75;2.3925;3.616666667];


[optp_ngpx,Popt_ngpx,~,ckdstpx] = findp_cont_ng(mc_px_ng,inp,Own);

sp_ng = optp_ngpx;
mc_ng = mc_px_ng;
pctmarkup_ng = (sp_ng-mc_ng)./sp_ng;
%%
inp.x2rest = rest_ng;
inp.dt3 = XS_ng2;
[optp_ngd,Poptwd3_ngd,dPcwd3_ngd] = findp_cont(mc_px_ng,inp,ownership);
pctmarkup_cf = (optp_ngd-mc_px_ng)./optp_ngd;

% [optp_ngpx_m2,Popt_ngpx_m2,~] = findp_cont_ng(mc_px_ng,inp,ownm2);
inp.x2rest = x2rest;
inp.dt3 = dt3;

%% Calulate utility change, wd2 
inp.ponly =0;
inp.inter1=0;
inp.wd2 = 1;
inp.wd3 = 0;

% U0 is observed utility before change
[~,~,U0] = predicted_share(thetafin_wd2,deltafin_wd2,inp);

Ub0 = U0.*JN;
ckub0 = zeros(1,R);
for i = 1:R
    i0=Ub0(:,i)==0;
    ckub0(1,i) = sum(i0);
end
sum(ckub0(:)==9);
ckckub0=R-sum(sum(JN));

Ub = sum(Ub0,1);

Ua0 = U_ngs.*PS_ngs;
Ua = sum(Ua0,1);
mui = -theta2wd2(2);

CSp = (Ub-Ua)./mui;

mCS = mean(CSp);

all(CSp(:)<0)

%% 
pchosel0=JN([1,4,7],:);
pchosel = sum(pchosel0,1);
ipl=pchosel';
age_l = zeros(4,2);
for d1 =1:4
    age_l(d1,1)=d1;
    ipl0 = age ==d1;
    age_l(d1,2)=sum(ipl0);
end

crosst_l = zeros(4,2);
for d1 =1:4
    ipl0 = age ==d1;
    djob0 = djob ==0;
    djob1 = djob ==1;
    ipl1 = ipl0.*djob0;
    ipl2 = ipl0.*djob1;
    crosst_l(d1,1)=sum(ipl1);
    crosst_l(d1,2)=sum(ipl2);
end
