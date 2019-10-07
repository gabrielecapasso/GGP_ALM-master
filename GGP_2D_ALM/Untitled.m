clear all
close all
load initial_random
nelx = 100;
MMC = true;
if MMC==true
        nely = 100;
        np = 5;                 % Number of deformable elements along x
        nY = 18;   % Number of deformable elements along y
        BC='L-shape';%L-shape %Short_Cantiliever
        p.method='MMC';%MMC%MNA %GP
        q=2;%q=1
        p.zp=1 ;% parameter for p-norm/mean regularization
        p.alp=1; %parameter for MMC
        p.epsi=0.7;% parameter for MMC
        p.bet=1e-3; %parameter for MMC
        p.deltamin=1e-6; %parameter for GP
        p.r=1.5;%parameter for GP
        minh=2;
        p.sigma=1.5;%parameter for MNA
        p.gammav=1;%parameter for GP
        p.gammac=3;%parameter for GP
        p.penalty=3;%parameter for MNA
        p.aggregation='KS'; %parameter for the aggregation function to be used
        % IE= Induced Exponential % KS= KS function %KSl= lowerbound KS function
        % p-norm %p-mean
        p.ka=4; % parameter for the aggregation constant
        p.saturation=false; % switch for saturation
        ncx=1; % number of components in the x direction
        ncy=1; % number of components in the y direction
        Ngp=1; % number of Gauss point per sampling window
        R=sqrt(3)/2; % radius of the sampling window (infty norm)
        initial_d=1; % component initial mass
end
Xc = Xg;
[W0,dW_dX,dW_dL,dW_dh]=Wgp_modified(ugp(:,1),ugp(:,2),Xc,p,np,nY,Yc,nely);
W = W0;
figure
imagesc(reshape(Aggregation_Pi(W,p),nelx,nely));
delta0=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
ddelta_dX=sum(reshape(full(dW_dX(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
ddelta_dL=sum(reshape(full(dW_dL(:,idgp).*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
ddelta_dh=sum(reshape(dW_dh(:,idgp).*repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3);

delta_c0=sum(reshape(W(:,idgp).^q.*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
if MMC==true
    ddelta_c_dX=sum(reshape(full(q*dW_dX(:,idgp).*repmat(W(:,idgp)'.^(q-1),1,nY*size(h,1))'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dL=sum(reshape(full(q*dW_dL(:,idgp).*repmat(W(:,idgp)'.^(q-1),1,nY*size(h,1))'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dh=sum(reshape(q*dW_dh(:,idgp).*repmat(W(:,idgp),size(h,1),1).^(q-1).*repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(h,1),1),size(h,1),[],Ngp^2),3);
else
    ddelta_c_dX=sum(reshape(full(q*dW_dX(:,idgp).*reshape(repmat(W(:,idgp)'.^(q-1),nY,1),size(W,2),[])'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dL=sum(reshape(full(q*dW_dL(:,idgp).*reshape(repmat(W(:,idgp)'.^(q-1),nY,1),size(W,2),[])'.*repmat(gauss_weight(:)',size(dW_dX,1),1)),size(dW_dX,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(dW_dX,1),1),size(dW_dX,1),[],Ngp^2),3);
    ddelta_c_dh=sum(reshape(q*dW_dh(:,idgp).*W(:,idgp).^(q-1).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
end
err = zeros(length(Xc)-np,1);
errd = err;
errdc = err;
var = 1e-3;
for i=1:2:2*np*nY
Xc1=Xc;
Xc1(i)=(1+var)*Xc(i);
dx = Xc1-Xc;
[W,~,~,~]=Wgp_modified(ugp(:,1),ugp(:,2),Xc1,p,np,nY,Yc,nely);
delta=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
delta_c=sum(reshape(W(:,idgp).^q.*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
err(i) = norm(W-W0-dx(1:2:2*np*nY)'*dW_dX)/norm(W0);
errd(i) = norm(delta-delta0-dx(1:2:2*np*nY)'*ddelta_dX)/norm(delta0);
errdc(i) = norm(delta_c-delta_c0-dx(1:2:2*np*nY)'*ddelta_c_dX)/norm(delta_c0);
end


for i=2:2:2*np*nY
Xc1=Xc;
Xc1(i)=(1+var)*Xc(i);
dx = Xc1-Xc;
[W,~,~,~]=Wgp_modified(ugp(:,1),ugp(:,2),Xc1,p,np,nY,Yc,nely);
delta=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
delta_c=sum(reshape(W(:,idgp).^q.*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
err(i) = norm(W-W0-dx(2:2:2*np*nY)'*dW_dL)/norm(W0);
errd(i) = norm(delta-delta0-dx(2:2:2*np*nY)'*ddelta_dL)/norm(delta0);
errdc(i) = norm(delta_c-delta_c0-dx(2:2:2*np*nY)'*ddelta_c_dX)/norm(delta_c0);
end
% 
% 
for i=2*np*nY+1:length(Xc)-np
Xc1=Xc;
Xc1(i)=(1+var)*Xc(i);
dx = Xc1-Xc;
[W,~,~,~]=Wgp_modified(ugp(:,1),ugp(:,2),Xc1,p,np,nY,Yc,nely);
delta=sum(reshape(W(:,idgp).*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
delta_c=sum(reshape(W(:,idgp).^q.*repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3)./sum(reshape(repmat(gauss_weight(:)',size(W,1),1),size(W,1),[],Ngp^2),3);
err(i) = norm(W-W0-dx(2*np*nY+1:length(Xc)-np)'*dW_dh)/norm(W0);
errd(i) = norm(delta-delta0-dx(2*np*nY+1:length(Xc)-np)'*ddelta_dh)/norm(delta0);
errdc(i) = norm(delta_c-delta_c0-dx(2*np*nY+1:length(Xc)-np)'*ddelta_c_dh)/norm(delta_c0);
% % norm(W-W0)
end
plot(err(err>0));
figure
plot(errd(errd>0));
figure
plot(errdc(errdc>0));