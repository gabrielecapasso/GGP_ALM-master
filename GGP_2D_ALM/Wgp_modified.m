function [W,dW_dX,dW_dL,dW_dh]=Wgp_modified(x,y,Xc,p,np,nY,Yk,nely)

Xk=Xc(1:2:2*np*nY);        % Extraire tous les Xk
Lk=Xc(2:2:2*np*nY);        % Extraire tous les Lk
h=Xc((2*np*nY)+1:end-np);
Xk=reshape(Xk,nY,np);
Lk=reshape(Lk,nY,np);
Yk=reshape(Yk,nY,np);
idy=ceil(y/max(Yk(:))*(nY-1));
idx=ceil(x/max(Xk(:))*(np-1));
% xk & dxk_dXk
xk=Xk(1:end-1,:);
dxk_dXk=zeros(size(xk,1),size(Xk,1));
dxk_dXk(1:end,1:(end-1))=eye(size(xk,1));
dxk_dXk=sparse(dxk_dXk);
% xk+1 & dxk+1_dXk
xk1=Xk(2:end,:);
dxk1_dXk=zeros(size(xk1,1),size(Xk,1));
dxk1_dXk(1:end,2:(end))=eye(size(xk1,1));
dxk1_dXk=sparse(dxk1_dXk);
% yk & yk+1
yk=Yk(1:end-1,:);
yk1=Yk(2:end,:);
% lk & dlk_dLk
lk=Lk(1:end-1,:);
dlk_dLk=zeros(size(lk,1),size(Lk,1));
dlk_dLk(1:end,1:(end-1))=eye(size(lk,1));
dlk_dLk=sparse(dlk_dLk);
% 1/lk & 1/dlk_dLk
d11lk_dLk=zeros(size(lk,1),size(Lk,1)); %%%
d11lk_dLk(1:end,1:(end-1))=-(1./(lk(1:end,1))).^2.*eye(size(lk,1));
d11lk_dLk=sparse(d11lk_dLk);
% lk+1 & dlk+1_dLk
lk1=Lk(2:end,:);
dlk1_dLk=zeros(size(lk,1),size(Lk,1));
dlk1_dLk(1:end,2:(end))=eye(size(lk1,1));
dlk1_dLk=sparse(dlk1_dLk);
switch p.method
    case 'MMC'
        ii=1:numel(x);
        jj=1:np;
        [I,J]=meshgrid(ii,jj);
        xg=reshape(x(I),size(I))';
        yg=reshape(y(I),size(I))';
        alp=p.alp;
        epsi=p.epsi;
        bet=p.bet;
        % ly & uy
        ly=xk(idy,:)-lk(idy,:)/2+(repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)).*((xk1(idy,:)-lk1(idy,:)/2)-(xk(idy,:)-lk(idy,:)/2));
        uy=xk(idy,:)+lk(idy,:)/2+(repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)).*((xk1(idy,:)+lk1(idy,:)/2)-(xk(idy,:)+lk(idy,:)/2));
        % derivé ly et uy 2704x90
        dly_Xk= repmat(dxk_dXk(idy,:),1,np)+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat((dxk1_dXk(idy,:)-(dxk_dXk(idy,:))),1,np);
        dly_Lk= repmat(-dlk_dLk(idy,:),1,np)/2+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat(((-dlk1_dLk(idy,:)/2)+dlk_dLk(idy,:)/2),1,np);
        duy_Xk= repmat(dxk_dXk(idy,:),1,np)+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat((dxk1_dXk(idy,:)-(dxk_dXk(idy,:))),1,np);
        duy_Lk= repmat(dlk_dLk(idy,:)/2,1,np)+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat((dlk1_dLk(idy,:)/2)-dlk_dLk(idy,:)/2,1,np);
        b=nely*h;
        b=repmat(b(:),1,size(y,1))';
        db_h =nely*ones(size(h,1),1); % 5x1
        xg = repmat(x,1,np);
        yg = repmat(y,1,np);
        % local chi
        chi0_fun = @(x,y) 1-(4*(x-(uy+ly)/2).^2./(uy-ly).^2).^alp-((y-b/2).^2./(1.05*b.^2)).^(alp);
        chi0 = chi0_fun(xg,yg)';
        dchi0_uy = (-(4*alp*(ly - xg).*((ly + uy - 2*xg).^2./(ly - uy).^2).^(alp - 1).*(ly + uy - 2*xg))./(ly - uy).^3)';
        dchi0_ly = ((4*alp*(uy - xg).*((ly + uy - 2*xg).^2./(ly - uy).^2).^(alp - 1).*(ly + uy - 2*xg))./(ly - uy).^3)';
        dchi0_b = (-alp*((20*(b/2 - yg))./(21*b.^2) - (40*(b/2 - yg).^2)./(21*b.^3)).*((20*(b/2 - yg).^2)./(21*b.^2)).^(alp - 1))';
%         dchi0_xg =(4*alp*((ly + uy - 2*xg).^2./(ly - uy).^2).^(alp - 1).*(ly + uy - 2*xg))./(ly - uy).^2;
%         dchi0_yg =(20*alp*(b - 2*yg).*((20*(b/2 - yg).^2)./(21*b.^2)).^(alp - 1))./(21*b.^2);
        [chi,dchi]=Aggregation_Pi(chi0,p);
        chi(chi<=-1e6)=-1e6; %1x2704
        dchi_uy = dchi0_uy.*dchi;
        dchi_ly = dchi0_ly.*dchi;
        dchi_b = dchi0_b.*dchi;
        W=(chi>epsi)+(chi<=epsi&chi>=-epsi).*(3/4*(1-bet)*(chi/epsi-chi.^3/3/epsi^3)+(1+bet)/2)+(chi<-epsi)*bet;
        dW_dchi=-3/4*(1/epsi-chi.^2/epsi^3).*(bet-1).*(abs(chi)<epsi);
        dW_chi=spdiags(dW_dchi',0,length(chi),length(chi));
        dW_chi = repmat(dW_dchi,nY*np,1);
        dW_chi2 = repmat(dW_dchi,np,1);
        dW_uy=dchi_uy.*dW_chi2;
        dW_ly=dchi_ly.*dW_chi2;
        dW_b=dchi_b.*dW_chi2;
        dW_dX=kron(dW_uy,ones(nY,1)).*duy_Xk'+kron(dW_ly,ones(nY,1)).*dly_Xk';
        dW_dL=kron(dW_uy,ones(nY,1)).*duy_Lk'+kron(dW_ly,ones(nY,1)).*dly_Lk';
        dW_dX=reshape(repmat(dW_uy,nY,1),length(idy),[])'.*duy_Xk'+reshape(repmat(dW_ly,nY,1),length(idy),[])'.*dly_Xk';
        dW_dL=reshape(repmat(dW_uy,nY,1),length(idy),[])'.*duy_Lk'+reshape(repmat(dW_ly,nY,1),length(idy),[])'.*dly_Lk';
        dW_dh=(repmat(db_h(:),1,size(y,1))).*dW_b;

    case 'GP'
        deltamin=p.deltamin;
        r=p.r;
        d1= xk(idy,:)- (repmat(y,1,np)+ lk(idy,:)/2);
        d2= -xk(idy,:)+ (repmat(y,1,np)- lk(idy,:)/2);
        q=10;
        zetavar= (1/q*log(0.5*(exp(q*d1) + exp(q*d2))));
        dzetavar_Xk_numinateur= 0.5*(repmat(exp(q*d1'),nY,1).*(q*repmat(dxk_dXk(idy,:),1,np))'- repmat(exp(q*d2'),nY,1).*(q*repmat(dxk_dXk(idy,:),1,np))');
        dzetavar_Xk_denominateur= 0.5*(repmat(exp(q*d1') + exp(q*d2'),nY,1)); 
        dzetavar_Xk= 1/q*(dzetavar_Xk_numinateur./dzetavar_Xk_denominateur);
        dzetavar_Lk_numinateur= 0.5*(repmat(exp(q*d1'),nY,1).*(q*repmat(-dlk_dLk(idy,:)/2,1,np))'- repmat(exp(q*d2'),nY,1).*(q*repmat(-dlk_dLk(idy,:)/2,1,np))');
        dzetavar_Lk_denominateur= 0.5*(repmat(exp(q*d1') + exp(q*d2'),nY,1)); 
        dzetavar_Lk= 1/q*(dzetavar_Lk_numinateur./dzetavar_Lk_denominateur);
        
        deltaiel=(1/pi/r^2*(r^2*acos(zetavar/r)-zetavar.*sqrt(r^2-zetavar.^2))).*(abs(zetavar)<=r)+((zetavar<-r));
        ddetlaiel_dzetavar=(-2*sqrt(r^2-zetavar.^2)/pi/r^2).*(abs(zetavar)<=r);
        W=deltamin+(1-deltamin)*deltaiel;
        dW_ddeltaiel=(1-deltamin);
        dW_dX=dW_ddeltaiel.* repmat(ddetlaiel_dzetavar',nY,1).*dzetavar_Xk;
        dW_dL=dW_ddeltaiel.*repmat(ddetlaiel_dzetavar',nY,1).*dzetavar_Lk;
        
    case 'MNA'
        % ly & uy
        ly=xk(idy,:)-lk(idy,:)/2+(repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)).*((xk1(idy,:)-lk1(idy,:)/2)-(xk(idy,:)-lk(idy,:)/2));
        uy=xk(idy,:)+lk(idy,:)/2+(repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)).*((xk1(idy,:)+lk1(idy,:)/2)-(xk(idy,:)+lk(idy,:)/2));
        % derivé ly et uy
        dly_Xk= repmat(dxk_dXk(idy,:),1,np)+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat((dxk1_dXk(idy,:)-(dxk_dXk(idy,:))),1,np);
        dly_Lk= repmat(-dlk_dLk(idy,:),1,np)/2+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat(((-dlk1_dLk(idy,:)/2)+dlk_dLk(idy,:)/2),1,np);
        duy_Xk= repmat(dxk_dXk(idy,:),1,np)+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat((dxk1_dXk(idy,:)-(dxk_dXk(idy,:))),1,np);
        duy_Lk= repmat(dlk_dLk(idy,:)/2,1,np)+reshape(repmat((repmat(y,1,np)-yk(idy,:))./(yk1(idy,:)-yk(idy,:)),nY,1),length(idy),[]).*repmat((dlk1_dLk(idy,:)/2)-dlk_dLk(idy,:)/2,1,np);
        gt=3;
        w =@(x) (0.5-(15/(16*gt))*x+(5/(8*gt^3))*x.^3-(3/(16*gt^5))*x.^5).*(x>=-gt&x<=gt)+(x<-gt);
        dw =@(x) (-15/(16*gt) +3*(5/(8*gt^3))*x.^2 -5*(3/(16*gt^5))*x.^4).*(x>=-gt&x<=gt);
        tv1=-repmat(x,1,np)+ly;
        tv2=repmat(x,1,np)-uy;
        b=nely*h;
        tv3=repmat(y,1,np)-(repmat(b(:),1,size(y,1)))';
        W=w(tv1)'.*w(tv2)'.*w(tv3)';
        %derivé tv1 et tv2
        dtv1_Xk=+dly_Xk;
        dtv2_Xk=-duy_Xk;
        dtv1_Lk=+dly_Lk;
        dtv2_Lk=-duy_Lk;
        %
        dh_h =ones(size(h,1),1);
        db_h=nely*ones(size(h,1),1);
        dtv3_h=-(repmat(db_h(:),1,size(y,1)))';
        %
        dW_dX =((reshape(repmat(dw(tv1).*w(tv2),nY,1),length(idy),[]).*dtv1_Xk+ (reshape(repmat(dw(tv2).*w(tv1),nY,1),length(idy),[]).*dtv2_Xk)).*(reshape(repmat(w(tv3),nY,1),length(idy),[])))' ;
        dW_dL =((reshape(repmat(dw(tv1).*w(tv2),nY,1),length(idy),[]).*dtv1_Lk+ reshape(repmat(dw(tv2).*w(tv1),nY,1),length(idy),[]).*dtv2_Lk) .*(reshape(repmat(w(tv3),nY,1),length(idy),[])))';
        dW_dh = ((reshape(dw(tv3),length(idy),[]).*dtv3_h).*(reshape((w(tv1).*w(tv2)),length(idy),[])))';
end