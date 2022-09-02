function [lklhdfn] = likelihd(theta,delta,inp)

%this function will be changed to likelihood function
%with given theta, it calculates the mainly the likelihood function, the case will be
%divided into two : concentrated or not concentrated
%for concentration case, another function optimization is needed.

JN=inp.JN ;


% calculate choice probabilities after concentration
[~,PS1] = predicted_share(theta,delta,inp);

logPS1 = log(PS1);
A = logPS1.*JN;

lklhdfn = -sum(sum(A));


end