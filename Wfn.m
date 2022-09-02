function [Wfn] = Wfn(inp,theta,Z,psdiff_fn)
  
J = inp.J;
R = inp.R;
u = psdiff_fn(theta);

if inp.pinter % When price in interaction terms, I should 
  ninst = size(Z,2);
  Wfn = zeros(ninst,ninst);
  for k = 1:R
      for j = 1:J
          Zu = Z(j,:)'*u(j,k);
          Wfn = Wfn +Zu*Zu';
      end
  end
   Wfn =Wfn/(J*R);
   rank(Wfn)
else 
   usqrd=u.*u;
   Wfn=sum(sum(usqrd));
   Wfn=Wfn/(J*R);
   rank(Wfn)
end
