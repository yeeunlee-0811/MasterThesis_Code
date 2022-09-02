function [fntominimize] = first_stage_obj(inp,Pdelta) 
 S = inp.S; % j by 1
 J = inp.J ;%the number of products (the number of possible combinations)
 fntominimize = 0;
 for k =1:J
     fntominimize = fntominimize + (Pdelta(k)-S(k))^2;
 end
end
     
 


