function beta=The_bayesian_lasso(x,y,n,p)




%betac_draw=zeros(6000,p);
%m=zeros(p,1);
%betac=mvnrnd(m,eye(p));
%betac=betac';
tau2=ones(p,1);
sigma2=0.1;


lambda=1;
%lambda=ones(p,1);
%lambda(seq)=0.01;
for iter=1:6000

%     for i=1:p
%         x_g=x(:,i);
%         x_dg=x;
%         x_dg(:,i)=[];
%         betac_g=betac(i);
%         betac_dg=betac;
%         betac_dg(i)=[];
%         sigmag=inv(x_g'*x_g+1/tau2(i));
%         miug=sigmag*x_g'*(y-x_dg*betac_dg);
% 
% 
%             betac(i)=normrnd(miug,sigmag*sigma2);
% 
%     end
    D=diag(tau2);
 
   sigmag =inv( x'*x + inv(D));

	betac = mvnrnd(sigmag*x'*y,sigma2*sigmag)';
   
  %  for i=1:p
 
%             pmu=lambda*sqrt(sigma2)/norm(betac(i));
%             plambda=lambda^2;
%             pd = makedist('InverseGaussian','mu',pmu,'lambda',plambda);
%             temp = random(pd);
%             if temp==0
%                 temp=0.01;
%             end
%             tau2(i)=1/temp;
%           end
        for i=1:p
par1=lambda^2;
par2=betac(i)^2/sigma2;
tau2(i)=gigrnd(0.5,par1,par2,1);

        end
            
      

    
    par1=p+0.01;
    par2=sum(tau2)/2+0.01;
    temp=sqrt(gamrnd(par1,1/par2));
%    lambda=ones(p,1)*temp;
  
  lambda=temp;
    
 %   lambda(seq)=0.01;
    
    D=diag(tau2);
%    temp=sum(gamma);
    p_alpha=n/2+p/2+0.01;
    p_beta=((y-x*betac)'*(y-x*betac)+betac'*inv(D)*betac)/2+0.01;
    temp2=gamrnd(p_alpha,1/p_beta);
    sigma2=1/temp2;
    


    beta_draw(iter,:)=betac;

end


beta=mean(beta_draw(2001:6000,:));



