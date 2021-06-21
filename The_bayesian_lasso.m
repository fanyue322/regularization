function z=Bayesian_lasso(x,y)
burn_in=3000;
burn_out=3000;
z=starty(0,1,y);
[n,p]=size(x);
tau2=ones(p,1);
sigma2=1;


m=zeros(p,1);
betac=mvnrnd(m,eye(p));

lambda=1;

for iter=1:3000

      D=diag(tau2);

   sigmag =inv( x'*x + inv(D));

	betac = mvnrnd(sigmag*x'*z,sigma2*sigmag)';
    
      for i=1:p
par1=lambda^2;
par2=betac(i)^2/sigma2;
tau2(i)=gigrnd(0.5,par1,par2,1);

      end
          par1=p+0.01;
    par2=sum(tau2)/2+0.01;
    temp=sqrt(gamrnd(par1,1/par2));
      lambda=temp;
      

    
          miu=x*betac;
    sigmaG=eye(n)*sigma2;
    z=multruncn(miu,sigmaG,z,y);
    
    
    
 beta_draw(iter,:)=betac;


end
z=mean(beta_draw(burn_in+1:burn_in+burn_out,:));

end
