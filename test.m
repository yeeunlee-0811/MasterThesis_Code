%t1 = {'1','2','3'}

%temp = str2double(t1)

%t2 =[13,4,5]

%t2(:,temp==1)=[]

%t3 = randn(3,3)
%t4 = randn(3,3)
%t5 = t4*inv(t3)
%t6 = t4*(t3\eye(3,3))

%t7 = eye(5,5)
%t8 = randn(5,1)
%t9 = t8.*t7
%t10 = sum(t9,1)'

%t11 = repmat(1:5,5,1)
%t12 = t11.*t10

%s1 =3
%t13 = repmat(t7(:,1),1,3)
%t14 =randn(1,5)
%t15 = t13.*t14

%t16 = reshape(1:25,5,5)
%t17 = [1,2,3,4,5]

%t17 = t14.*repmat(t16(:,1),1,5)

%t3 = randn(3,3)
%t18 = t3.*t3'
%all(eig(t18)>eps)
%[~,p] = chol(t18)

%t16 = reshape(1:25,5,5)
%t20 =diag(t16)
%t19 = t16'*t16

%t21 = t3\eye(3,3)
%t22 =inv(t3)
%t23 = zeros(1,size(row,1));
%for i =1:size(row,1)
%find(schoice ==row(i,1))
%end 

t24 = randn(5,5)
% t25 = randn(5,1)
% 
% t26=t24./t25
% 
% Ipre = [1;1;1;2;2;2;3;3;3];
% I = zeros(J,J)
% for i =1:J
%     I(:,i) = Ipre ==Ipre(i)
% end

t25 = t24+5
t26=1./t25

t27 = eye(5,5)
t24.*t27(1,:)