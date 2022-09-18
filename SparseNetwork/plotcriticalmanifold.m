% clear all; 

Delta0=0.3;J0=1.0;taum=1;M=1e+3;
etabar=-0.5; 
DeltaEta=1e-4;
%%% alpha=-1.0 : inhibitory case (diVolo/Torcini)
%%% alpha= 1.0 : excitatory case
alpha=1.0; 


V = linspace(-6,-0.1,100000);
R=-DeltaEta./((2*V*pi+J0*Delta0)*taum);
K=(-V.^2+pi^2*taum^2*R.^2)/sqrt(M)-alpha*J0*taum*R;
plot3(K,V,R,'-','color',[0.7 0.7 0.7],'linewidth',3);hold on;

V = linspace(-0.1,-0.001,1000000);
R=-DeltaEta./((2*V*pi+J0*Delta0)*taum);
K=(-V.^2+pi^2*taum^2*R.^2)/sqrt(M)-alpha*J0*taum*R;
plot3(K,V,R,'-','color',[0.7 0.7 0.7],'linewidth',3);hold on;

V = linspace(0.001,0.1,1000000);
R=-DeltaEta./((2*V*pi+J0*Delta0)*taum);
K=(-V.^2+pi^2*taum^2*R.^2)/sqrt(M)-alpha*J0*taum*R;
plot3(K,V,R,'-','color',[0 0 1],'linewidth',3);hold on;

% V = linspace(0.1,6,1000000);
% R=-DeltaEta./((2*V*pi+J0*Delta0)*taum);
% K=(-V.^2+pi^2*taum^2*R.^2)/sqrt(M)-alpha*J0*taum*R;
% plot3(K,V,R,'-','color',[0 1 0],'linewidth',3);hold on;
% hold off;

% hold on;
% interval=[-1 1 -6 8 0 200];
% % g = @(K,V,R) DeltaEta/(pi*taum) + 2*R.*V + R*J0*Delta0/pi;
% g = @(K,V,R) V.^2 + sqrt(M)*(K+J0*taum*R)-(pi*taum*R).^2;
% fimplicit3(g,interval,'FaceColor',[0.7 0.7 0.7],'EdgeColor','none','FaceAlpha',.5);
% hold off;

% hold on;
% interval=[-1 1 -6 8 0 200];
% vStar = -J0*Delta0/(2*pi);
% g = @(K,V,R) vStar.^2 + sqrt(M)*(K+J0*taum*R)-(pi*taum*R).^2;
% fimplicit3(g,interval,'FaceColor',[1.0 0.0 0.0],'EdgeColor','none','FaceAlpha',.5);
% hold off;

% axis(interval);



%%%
% set(gcf,'color','white');
% set(gca,'FontName','Times','FontSize',14);
% axis([-1 1 -1 1 -8 8]);
% set(gca,'XTick',-1:1:1);
% set(gca,'YTick',-1:1:1);
% set(gca,'ZTick',-8:8:8);
% view([22 13]);
% set(gcf,'color','white');
 %set(gca,'FontName','Times','FontSize',14);
% axis([-1 1 -20 20 -6 6]);
% set(gca,'XTick',-1:1:1);
% set(gca,'YTick',-20:20:20);
% set(gca,'ZTick',-6:6:6);

