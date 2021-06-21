

function z=bayesian_fused_lasso(x,y,n,p)
z=starty(0,1,y);
[n,p]=size(x);

     burn_in=3000;
    burn_out=3000;
beta_c=normrnd(0,1,p,1);
tau2=ones(p,1);
for i=1:p
    tau2(i)=gamrnd(1,1);
end
ome2=ones(p-1,1);
for i=1:p-1
   ome2(i)=gamrnd(1,1);
 % ome2(i)=inf;
end
lambda1=1;
lambda2=1;
sigma2=1;
beta_draw=zeros(burn_in+burn_out,p);


    for iter=1:burn_in+burn_out
        sigmab=zeros(p,p);
        for i=1:p
            if i>1&&i<p
                sigmab(i,i)=1/tau2(i)+1/ome2(i)+1/ome2(i-1);
            elseif i==1
                sigmab(i,i)=1/tau2(i)+1/ome2(i);
            elseif i==p
                sigmab(i,i)=1/tau2(i)+1/ome2(i-1);
            end
        end
        for i=1:p-1
            sigmab(i+1,i)=-1/ome2(i);
        end
        for i=1:p-1
            sigmab(i,i+1)=-1/ome2(i);
        end
       mug=inv(x'*x+sigmab)*x'*z;
%mug=(x'*x+sigmab)\x'*y;
     
      sigmag=sigma2*inv(x'*x+sigmab);
%      temp=sigma2*eye(p);
     % sigmag=temp*inv(x'*x+sigmab);
   %   sigmag=temp/(x'*x+sigmab);
%        sigmag=sigma2/(x'*x+sigmab);
        beta_c=mvnrnd(mug,sigmag);
        beta_c=beta_c';
        for i=1:p
%             pmu=sqrt(lambda1^2*sigma2/beta_c(i)^2);
%             plambda=lambda1^2;
%             pd = makedist('InverseGaussian','mu',pmu,'lambda',plambda);
%             temp = random(pd);
%             tau2(i)=1/temp;
par1=lambda1^2;
par2=beta_c(i)^2/sigma2;
tau2(i)=gigrnd(0.5,par1,par2,1);
        end
        
        for i=1:p-1
%             pmu=sqrt(lambda2^2*sigma2/(beta_c(i+1)-beta_c(i))^2);
%             plambda=lambda2^2;
%             pd = makedist('InverseGaussian','mu',pmu,'lambda',plambda);
%             temp = random(pd);
%             ome2(i)=1/temp;
par1=lambda2^2;
par2=(beta_c(i+1)-beta_c(i))^2/sigma2;
ome2(i)=gigrnd(0.5,par1,par2,1);
        end
        
        
        sigmab=zeros(p,p);
        for i=1:p
            if i>1&&i<p
                sigmab(i,i)=1/tau2(i)+1/ome2(i)+1/ome2(i-1);
            elseif i==1
                sigmab(i,i)=1/tau2(i)+1/ome2(i);
            elseif i==p
                sigmab(i,i)=1/tau2(i)+1/ome2(i-1);
            end
        end
        for i=1:p-1
            sigmab(i+1,i)=-1/ome2(i);
        end
        for i=1:p-1
            sigmab(i,i+1)=-1/ome2(i);
        end
        
        
        
        
        
        
        par1=(n+p-1)/2+0.01;
        par2=0.5*(z-x*beta_c)'*(z-x*beta_c)+0.5*beta_c'*sigmab*beta_c+0.01;
       temp=gamrnd(par1,1/par2);
       sigma2=1/temp;
       
       par1=p+0.1;
       par2=0.5*sum(tau2)+0.1;
       lambda1=gamrnd(par1,1/par2);
       lambda1=sqrt(lambda1);
       
       par1=p-1+0.1;
       par2=0.5*sum(ome2)+0.1;
       lambda2=gamrnd(par1,1/par2);
       lambda2=sqrt(lambda2);
       
       
     miu=x*betac;
    sigmaG=eye(n)*sigma2;
    z=multruncn(miu,sigmaG,z,y);

       beta_draw(iter,:)=beta_c;
    end
    

    
    
    

z=mean(beta_draw(burn_in+1:burn_in+burn_out,:));



end

        
            
