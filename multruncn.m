function Y=multruncn(mu,S,Y,Z)

% multruncn: Bayesian Variable Selection
%       Sampling from multivariate truncated normal 
%       based on class label Z,
%       with mean vector mu and covariance matrix S 
%       Nominal response case.
%
%Copyright (c) 2003 Naijun Sha.
%**********************************************************************


n=size(Y,1);
for i=1:n
   Y2=Y;Y2(i)=[];
   mu1=mu(i);mu2=mu;mu2(i)=[];
   S11=S(i,i);S12=S(i,:);S12(:,i)=[];S21=S12';
   S22=S;S22(i,:)=[];S22(:,i)=[];S22inv=inv(S22);
   S112=S11-S12*S22inv*S21;s=sqrt(S112);
   m1=mu1+S12*S22inv*(Y2-mu2);p=normcdf(0,m1,s);u=rand(1);
   z=Z(i);
   switch z
   case 0
      Y1=norminv(u*p,m1,s);
   case 1
      Y1=norminv(u*(1-p)+p,m1,s);
   end;
   Y(i)=Y1;
end;



   








