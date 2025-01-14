function [almc,dalmc,mtv] = ALM_constraint(Xc,nely,nY,np,p)

Xk=Xc(1:2:end-2*np);       
Lk=Xc(2:2:end-2*np);        
Xk=reshape(Xk,nY,np);
Lk=reshape(Lk,nY,np);

j1=reshape(repmat(1:np,nY-1,1),[],1);
j2=nY*(repmat(j1,1,nY)-1)+repmat(1:nY,length(j1),1);
i2=repmat((1:(np*(nY-1)))',1,nY);
% xk & dxk_dXk
xk=Xk(1:end-1,:);
dxk_dXk=zeros(size(xk,1),size(Xk,1));
dxk_dXk(1:end,1:(end-1))=eye(size(xk,1));
dxk_dXk=repmat(dxk_dXk,np,1);
dxk_dXk=sparse(i2(:),j2(:),dxk_dXk(:),(nY-1)*np,nY*np);
% xk+1 & dxk+1_dXk
xk1=Xk(2:end,:);
dxk1_dXk=zeros(size(xk1,1),size(Xk,1));
dxk1_dXk(1:end,2:(end))=eye(size(xk1,1));
dxk1_dXk=repmat(dxk1_dXk,np,1);
dxk1_dXk=sparse(i2(:),j2(:),dxk1_dXk(:),(nY-1)*np,nY*np);
% lk & dlk_dLk
lk=Lk(1:end-1,:);
dlk_dLk=zeros(size(lk,1),size(Lk,1));
dlk_dLk(1:end,1:(end-1))=eye(size(lk,1));
dlk_dLk=repmat(dlk_dLk,np,1);
dlk_dLk=sparse(i2(:),j2(:),dlk_dLk(:),(nY-1)*np,nY*np);
% lk+1 & dlk+1_dLk
lk1=Lk(2:end,:);
dlk1_dLk=zeros(size(lk,1),size(Lk,1));
dlk1_dLk(1:end,2:(end))=eye(size(lk1,1));
dlk1_dLk=repmat(dlk1_dLk,np,1);
dlk1_dLk=sparse(i2(:),j2(:),dlk1_dLk(:),(nY-1)*np,nY*np);

uuk=xk+lk/2;
uuk1=xk1+lk1/2;
llk=xk-lk/2;
llk1=xk1-lk1/2;
duuk_dXk=dxk_dXk;
duuk_dLk=dlk_dLk/2;
duuk1_dXk=dxk1_dXk;
duuk1_dLk=dlk1_dLk/2;
dllk_dXk=dxk_dXk;
dllk_dLk=-dlk_dLk/2;
dllk1_dXk=dxk1_dXk;
dllk1_dLk=-dlk1_dLk/2;

tanalpha= ((uuk1 - uuk)/nely)*(nY-1);
tanbeta= ((llk - llk1)/nely)*(nY-1);
dtanalpha_dXk= ((duuk1_dXk - duuk_dXk)/nely)*(nY-1);
dtanbeta_dXk= ((dllk_dXk - dllk1_dXk)/nely)*(nY-1);
dtanalpha_dLk= ((duuk1_dLk - duuk_dLk)/nely)*(nY-1);
dtanbeta_dLk= ((dllk_dLk - dllk1_dLk)/nely)*(nY-1);
tv1 = ((tanalpha(:) - tand(45))/tand(45));
tv2= ((tanbeta(:) - tand(45))/tand(45));
TV=[tv1;tv2];
mtv=max(TV);
p.aggregation='KS';
[almc,dalmc_dtv]=Aggregation_Pi(TV,p);
almc=almc*100;

dtv1_dXk = ((dtanalpha_dXk)/tand(45));
dtv1_dLk = ((dtanalpha_dLk)/tand(45));
dtv2_dXk = ((dtanbeta_dXk)/tand(45));
dtv2_dLk = ((dtanbeta_dLk)/tand(45));

dtv_dXk=[dtv1_dXk;dtv2_dXk];
dtv_dLk=[dtv1_dLk;dtv2_dLk];
dalmc_dXk=dalmc_dtv(:)'*dtv_dXk;
dalmc_dLk=dalmc_dtv(:)'*dtv_dLk;
dalmc=zeros(size(Xc));
dalmc(1:2:(end-2*np))=100*dalmc_dXk;
dalmc(2:2:(end-2*np))=100*dalmc_dLk;

end

