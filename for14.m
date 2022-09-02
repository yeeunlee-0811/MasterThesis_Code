%Drop individuals : 2423 Obs.
% XIinter_indiv = XIinter;
% XIinter_indiv(row_indiv,:) =[];
% 
% inp.XI_indiv = XIinter_indiv;
% inp.Rindiv = Rindiv;
% 
% % JN_indiv
% redschoice_indiv = redschoice;
% Wnnchose = redschoice(row_indiv);
% Sadjust = zeros(9,1);
% for k = 1:9
%     Sadjust(k,1) = sum(Wnnchose(:)==k);
% end
% Snumobs = S.*R;
% Snumindiv = Snumobs - Sadjust;
% Sindiv = Snumindiv./Rindiv;
% logSindiv = log(Sindiv);
% inp.logSindiv = logSindiv;
% inp.Sinidv= Sindiv;
% 
% redschoice_indiv(row_indiv,:) = [];
% 
% JN_indiv = zeros(J,Rindiv);
%  for j = 1:J
%      for r=1:Rindiv
%          if redschoice_indiv(r,1)==j
%          JN_indiv(j,r) = 1;
%          end
%      end
%  end
%  
% inp.JN_indiv = JN_indiv;
% 
% Xindiv1 = repmat(Xindiv1,1,J);
% inp.Xindiv1 = Xindiv1;
% 
% Xindiv2 = repmat(Xindiv2,1,J);
% inp.Xindiv2 = Xindiv2;
% 
% nindiv = 3*2;
% inp.nindiv = nindiv;