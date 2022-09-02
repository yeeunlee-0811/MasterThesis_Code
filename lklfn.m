function [lklhoodfn] = lklfn(P,inp)
%% input description
% indicator matrix
JN =inp.JN;
% the number of individual
R = inp.R;
% the number of Products
J = inp.J;
%% generating likelihood function
        lklhoodfn = -prod((P.^JN),'all');
 
end
    