clear
%% Importing data
% Importing individual data
panelm = importdata('exp_est1.xls');
XId = panelm;
indiv_varnames = XId(1,:);
XId = array2table(XId, 'VariableNames',indiv_varnames);
XId(:,1)=[]; % drop person ID
XId(1,:)=[]; % drop varnames

% Importing Service data
XSd = importdata('spconv.xls');
XS = XSd.data;
XSvarnames = XSd.textdata(1,:) ;
XSvarnames(:,1)=[];
XSd = array2table(XS, 'VariableNames', XSvarnames);

% Importing Phone data
%XP = importdata('pp.xls');
%ptype = XP.textdata(:,1);
%ptype(1,:) =[];
%XPd=XP.data;
%XPd = [ptype num2cell(XPd)];
%XPvarnames = XP.textdata(1,:);
%XP = array2table(XPd, 'VariableNames', XPvarnames);

%% Data cleaning for estimation 
% Individual data

XIdd=XId(:,{'hhldsiz15','age115', 'school15', 'gender','djob'});
XI = table2array(XIdd);
XI = str2double(XI);

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

[R,nikrc]=size(XI);

% Service data
XS = table2array(XSd);
%%%% remove products whose shares are zero
row = find(XS(:,19) ==0);
XS(row,:) =[]; 

Sfirm = XS(:,20:22);
sfirmid=zeros(size(Sfirm,1),1);
for j = 1:3 
    sfirmid=sfirmid + Sfirm(:,j)*j;
end

S = XS(:,19);

XS2=XS;
XS2(:,[10:17,19:23])=[];
[J,nskrc] = size(XS2);

% jn matrix
JN = zeros(J,R);
 for j = 1:J
     for r=1:R
         if schoice(r,1)==j
         JN(j,r) = 1;
         end
     end
 end
chkJN =zeros(1,R);
 for i = 1:R
chkJN(i)=sum( JN(:,i)==1 );
 end
%all(chkJN ==1)

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
dage = zeros(size(XId,1),3);
for k=1:3
    ida = age == k+1;
    dage(:,k)=ida;
end
agenpay = dage.*npayer;

XIinter = [incpay,agenpay]; % R by 3+3

XSinter= repmat(XS(:,18),1,6); %  J by 3+3

ninter =size(XIinter,2);
if size(XIinter,2) ~=size(XSinter,2)
     disp('The number of interaction terms are not matched')
else
    ninter = size(XIinter,2);
end


logS = zeros(size(S,1),1);
for i = 1:size(S,1)
    if S(i) ~=0
        logS(i) = log(S(i));
    else 
        logS(i) = log(1e-100);
    end
end
    
%Xdp = eye(J,J);
%drop the base product : 3G
%Xdp(:,14)=zeros(J,1);
% As a result, Xdp is a J by J-1 MATRIX 
%% Input Structure

inp.R = R;
inp.nikrc = nikrc;
inp.J = J;
inp.nskrc=nskrc;
inp.JN = JN;
inp.S= S;
inp.ninter = ninter; 
inp.logS = logS;

inp.XI = XI'; %%%%
inp.XS2 = XS2;
inp.XIinter = XIinter;
inp.XSinter = XSinter;
%inp.Xdp = Xdp;
%% 
obj = @(theta) lklfn2(inp,theta);

options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-6);

theta0 = 0.5*ones(ninter,1);
thetafin = fminunc(obj, theta0, options);
%% deltafin
[~,~,~,deltafin] = lklfn2(inp,thetafin);
[Pfin] = predicted_share(thetafin,deltafin,inp);

ntheta = size(thetafin,1); % ninter + nikrc
parmnum1 = ntheta + J;
%% 2nd stage parameter 
%% Second stage under likelihood method :
% Instruments: Interaction terms & Service price (before discounts)
spbd = XS(:,10)/10000;

%changing the order of variables
ssp =XS2(:,10);
XS20=XS2;
XS20(:,10)=[];
XS2=[ones(J,1),ssp, XS20];
% Explanatory variables
X=XS2;
INS = [ones(J,1),spbd,XS20];

[theta1,stderr,pval]=tsls(deltafin, X,INS);
%% Inputr Structure

inp.delta = deltafin;
inp.theta = thetafin;
inp.share = Pfin;
inp.eps = 6.0554544523933429e-6;
inp.ntheta = ntheta;
inp.parmnum1 = parmnum1;
inp.Dlikelmat0 = zeros(R,ntheta);

[dh] = stepP(inp);
inp.dh = dh;

[dhd] = stepD(inp);
inp.dhd = dhd;

%% Calculating se
[V1, Dlikelmat]=varMOM(inp);
%% 
Gam = Gamma(inp, Dlikelmat);
%% 
%G1 = Gam(1:J,1:J);
%V3 = V1(1:J,1:J);
%VQ = (inv(G1'*G1)*G1'*V3*G1*inv(G1'*G1))./R;
%sed = diag(sqrt(VQ));

G2 = Gam(J+1:J+ntheta, J+1:J+ntheta);
V4 = V1(J+1:J+ntheta,J+1:J+ntheta);
VW = inv(G2'*G2)*G2'*V4*G2*inv(G2'*G2)./R;
sew = diag(sqrt(VW));
