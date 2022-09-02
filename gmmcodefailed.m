%% GMM Method

% To write GMM, I should make instruments first
% Demand side instrument: 
instr_D = XSd;
instr_D(:,[10:18,20:end])=[]; 
instr_D = table2array(instr_D);
row = find(instr_D(:,end) ==0); 
instr_D(row,:) =[]; % drop products whose shares are zero
instr_D(:,end)=[];
L_D = size(instr_D,2);

z_D = zeros(J,L_D,3);
for i=1:J
     same_firm = sfirmid==sfirmid(i) ;
     same_firm
     diff_firm = sfirmid~=sfirmid(i) ;
     diff_firm
    for l=1:L_D
        zk = instr_D(:,l);
        z_D(i,l,1) = zk(i);
        z_D(i,l,2) = sum(zk(same_firm)) - zk(i);  
        z_D(i,l,3) = sum(zk(diff_firm));
    end
end
z_D = reshape(z_D,J,L_D*3);

z0_D = z_D(:,1:2);
for i=3:size(z_D,2)
    zi = z_D(:,i);
    i
    reg = fitlm(z0_D,zi,'Intercept',false);
    if reg.Rsquared.Ordinary < 0.98
        if length(unique(zi))>1 % not constant
            zi = zi - mean(zi); %?? why subtract mean ??%
        end
        z0_D = [z0_D zi];
        % disp(i)
    end
end

z_D = z0_D;

inp.zD=z_D;
%% 
W = (z_D'*z_D)/J; % initial weighting matrix
inp.W = W;
inp.conc = true; % "concentrate out" linear parameters
%inp.supply_side = false; % do *not* estimate supply side 
%%
psdiff_fn = @(theta) ps_diff(inp,theta);
inp.instoptimal = true;
inp.pinter = false;
obj = @(beta) GMM_objective(beta,z_D,W,psdiff_fn,inp);
theta0 = 0.5*ones(ninter,1); % starting values; randomly drawn from (0,1)

options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
theta2_0 = fminunc(obj,theta0,options); % theta2 = (sigma,alpha)
% options = optimset('Display','iter','MaxFunEvals',20000);
% theta2_0 = fminsearch(obj,theta0,options);
%%
[~,beta_0] = residfun(inp,theta2_0); % get linear parameters, 
thetagmm = [beta_0; theta2_0];
%% 
heterosked = true;
linear = false;
inp.conc = false;
psdiff_fn = @(theta) ps_diff(inp, theta);
stdErr = GMM_se(theta2_0,z_D,inp,psdiff_fn,linear,heterosked);
disp('Demand estimates')
disp([theta2_0, stdErr])

%% Second stage under likelihood method :
% Instruments: Firm dummy (base : lg)----------------------------What else?
zsls = Sfirm;
zsls(:,end)=[];
% making explanatory vars
XSmp=XS;
P = XS(:,end);
XSmp(:,end)=[];
x1s = [zsls,XSmp];
rank(x1s)
% prediction
sreg1 = fitlm(x1s,P,'Intercept',false);
fittedP = sreg1.Fitted;
% 2ND stage
x2s = [XSmp,fittedP];
sreg2 = fitlm(x2s,deltafin,'Intercept',false)

%% <<<<<<<<<<<<<<<<<<<<<<<<< GMM method >>>>>>>>>>>>>>>>>>>>>>>>>>
%% Instruments : Interaction terms : data * age115  / BLP instruments ?
% Interaction term as instruments : Interaction term mean & firm dummy
Zd_inter = zeros(9:3);
Zd_inter1 = repmat(XIinter',J,1).*JN;
Zd_inter(:,1)=XSinter.*mean(Zd_inter1,2); %%%--------------------------------> Is this right?
Zd_inter(:,2) = [ones(3,1);zeros(6,1)];
Zd_inter(:,3)=[zeros(3,1);ones(3,1);zeros(3,1)];
% BLP Instruments 
sfirmid =ones(9,1);
sfirmid = sfirmid +[zeros(3,1);ones(6,1)];
sfirmid = sfirmid +[zeros(6,1);ones(3,1)];
Zd_BLP = zeros(J,3);
instr_D = XSred(:,1);
for i=1:J
     same_firm = sfirmid==sfirmid(i) ;
     diff_firm = sfirmid~=sfirmid(i) ;
       zk = instr_D;
        Zd_BLP(i,1) = zk(i);
        Zd_BLP(i,2) = sum(zk(same_firm)) - zk(i);  
        Zd_BLP(i,3) = sum(zk(diff_firm));
end

z_D = Zd_inter;
inp.zD = z_D;
%% GMM parameter estimation
W = (z_D'*z_D)/J; % initial weighting matrix
inp.W = W;
inp.conc = true; % "concentrate out" linear parameters
psdiff_fn = @(theta) ps_diff(inp,theta);
inp.pinter = false;
obj = @(beta) GMM_objective(beta,z_D,W,psdiff_fn,inp);
theta0 = 0.5*ones(ninter,1); % starting values; randomly drawn from (0,1)

options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-6);
theta2_0 = fminunc(obj,theta0,options); % theta2 = (sigma,alpha)
% options = optimset('Display','iter','MaxFunEvals',20000);
% theta2_0 = fminsearch(obj,theta0,options);
[~,beta_0] = residfun(inp,theta2_0); % get linear parameters, 
thetagmm = [beta_0; theta2_0];
%% GMM se estimation
W = (z_D'*z_D)/J; % initial weighting matrix
inp.W = W;
heterosked = true;
linear = false;
inp.conc = false;
psdiff_fn = @(theta) ps_diff(inp, theta);
se0gmm = GMM_se(thetagmm,z_D,inp,psdiff_fn,linear,heterosked);
disp('Demand estimates')
disp([thetagmm stdErr])