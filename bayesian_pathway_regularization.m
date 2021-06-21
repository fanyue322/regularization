function beta_draw=bayesian_pathway_regularization(x,y,L)

burn_in=1000;
burn_out=2000;
r=1;
z=starty(0,1,y);
[n,p]=size(x);

tau2=ones(p,1);
 beta=zeros(p,1);
 sigma2=1;
 lambda=1;

D=eye(p);

prA=r*(inv(D)+L{1});


for iter=1:burn_in+burn_out

  
    miub=inv(x'*x+prA)*x'*z;
    sigmab=sigma2*inv(x'*x+prA);
    beta=mvnrnd(miub,sigmab);
    beta=beta';
 
      

        for i=1:p

par1=lambda^2;
par2=r*beta(i)^2/sigma2;
tau2(i)=gigrnd(0.5,par1,par2,1);



        end
   D=diag(tau2);
prA=r*(inv(D)+L{1});

% par1=(n+p)/2+0.01;
% par2=((z-x*beta)'*(z-x*beta)+beta'*prA*beta)/2+0.01;
% temp=gamrnd(par1,1/par2);
% sigma2=1/temp;

par1=p+0.01;
par2=sum(tau2)+0.01;
lambda=gamrnd(par1,1/par2);
 lambda=sqrt(lambda);

par1=p/2+0.01;
par2=beta'*prA*beta/(2*sigma2)+0.01;
r=gamrnd(par1,1/par2);
r=1;
D=diag(tau2);
prA=r*(inv(D)+L{1});

          miu=x*beta;
    sigmaG=eye(n)*sigma2;
    z=multruncn(miu,sigmaG,z,y);

    beta_draw(iter,:)=beta;
end  


%beta= mean(beta_draw(burn_in+1:burn_in+burn_out,:));