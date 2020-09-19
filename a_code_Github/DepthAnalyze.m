function [BathyFeatures]=DepthAnalyze(East,North,Depth,DepthCenter)
% data assignment
xfit=East-(max(East)+min(East))*1.0/2;
yfit=North-(max(North)+min(North))*1.0/2;
zfit=Depth-0;
Zarea=Depth;
Zgrid=DepthCenter(:,3);
%fitting function (Evan's)
func = @(var,x) var(1)*x(:,1).^2 + var(2)*x(:,2).^2 + var(3)*x(:,1).*x(:,2) + var(4)*x(:,1) + var(5)*x(:,2) + var(6);
a = lsqcurvefit(func,ones(1, 6),[xfit, yfit],zfit);
% computing the fitting errors
Zfitting= a(1) * xfit.^2 + a(2) * yfit.^2 + a(3) * xfit.*yfit + a(4) * xfit + a(5) * yfit + a(6);
DetZ=zfit-Zfitting;
%% bathymetry data analysis
av=a(1);
bv=a(2);
cv=a(3);
dv=a(4)+0.00001;
ev=a(5)+0.00001;
fv=a(6);
% (1)slope
SlopeR=atan(sqrt(dv*dv+ev*ev)); % radian
SlopeA=SlopeR*1.0/pi*180;
% (2)aspect,eastness,northness
AspectD=atan(ev/dv);
EastD=sin(AspectD);
NorthD=cos(AspectD);
% (3)curvature,BPI
ProCur=-200.0*(av*dv*dv+bv*ev*ev+cv*dv*ev)/((ev*ev+dv*dv)*((1+ev*ev+dv*dv).^1.5));
PlanCur=200.0*(bv*dv*dv+av*ev*ev-cv*dv*ev)/((ev*ev+dv*dv).^1.5);
ProCmax=-av-bv+sqrt((av-bv).^2+cv*cv);
ProCmin=-av-bv-sqrt((av-bv).^2+cv*cv);
MeanCur=-av-bv;
BPI=Zgrid-mean(Zarea);
% (4)TRI, roughness
TRI=sum(abs(Zarea-Zgrid))/(length(Zarea)-1);
RoughD=max(Zarea)-min(Zarea);
% (5) rugosity (plane and E-N projected)
% xyz=[xfit yfit zfit];
% tri = delaunay(xfit,yfit);
% [RugosityPCA, s, a, ~, ~, ~, RugosityXY] = trisurfterrainfeats(tri, xyz);
%% Bathymetric features
% BathyFeatures=[SlopeR, EastD,NorthD, ProCur PlanCur ProCmax ProCmin MeanCur  BPI TRI RoughD];% RugosityPCA RugosityXY];
BathyFeatures=[BPI TRI RoughD ProCur PlanCur ProCmax ProCmin MeanCur SlopeR,AspectD, EastD, NorthD];
end
