clear

J=10;
R=100;
nskrc=20;
ninter=5; % The number of interaction terms, first five Krc of Krc with individual preference
nikrc = 15;
% Suppose there exists three instruments
Z= randn(J,3);

inp.J = J;
inp.R = R;
% the number of product characteristics
inp.nskrc= nskrc;
inp.ninter = ninter;
inp.XI = randn(nikrc,R) ;
inp.XS = randn(J,nskrc);
inp.JN = repmat(eye(10),1,10);
inp.S = rand(J,1);
inp.nikrc =nikrc;
inp.Z =Z;
inp.Xinter= randn(ninter,R);
obj = @(theta) likelihd(theta,[],inp);

options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',20000,...
    'OptimalityTolerance',1e-6,'StepTolerance',1e-6);

theta0 = 0.5*ones(ninter+nikrc,1);
thetafin = fminunc(obj, theta0, options);
[~,deltafin] = likelihd(thetafin,[],inp);
[Pfin] = predicted_share(thetafin,deltafin,inp);


