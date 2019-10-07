close all
clc
x = linspace(0,2.5,101);
[xx,yy] = meshgrid(x,20*x);
a=3;
z1=1-((xx-1).^2).^a-((yy-20).^2/410).^a;
z2 = (1-((xx-1).^2).^a).*(1-((yy-20).^2/400).^a);
bet=1e-3;
epsi=0.7;
chi = z1;
W=(chi>epsi)+(chi<=epsi&chi>=-epsi).*(3/4*(1-bet)*(chi/epsi-chi.^3/3/epsi^3)+(1+bet)/2)+(chi<-epsi)*bet;
chi = z2;
W2=(chi>epsi)+(chi<=epsi&chi>=-epsi).*(3/4*(1-bet)*(chi/epsi-chi.^3/3/epsi^3)+(1+bet)/2)+(chi<-epsi)*bet;
figure
surfc(xx,yy,W);
figure
surfc(xx,yy,W2);
figure
contour(xx,yy,W);
figure
contour(xx,yy,W2);
%%
syms xg uy ly yg b alp
chi = 1-(4*(xg-(uy+ly)/2).^2./(uy-ly).^2).^alp-((yg-b/2).^2./(1.05*b.^2)).^(alp);
dchi0_uy = simplify(diff(chi,uy))
dchi0_ly = simplify(diff(chi,ly))
dchi0_b = simplify(diff(chi,b))
dchi0_xg = simplify(diff(chi,xg))
dchi0_yg = simplify(diff(chi,yg))
