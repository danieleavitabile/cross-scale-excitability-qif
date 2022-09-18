clear all, close all;

n      = 1e4;           M = 1e3;
J      = 1/sqrt(M); tauM  = 1;
vt     = 100;        vr   = -vt; 
tauD   = 1.5e-2;   Delta0 = 0.3; 
eps = 0.1; %0.5;
etaBar   = -0.5*sqrt(M);
DeltaEta = 1e-4*sqrt(M);
idPerm   = randperm(n)';
saveFiles = true;
T = 2*pi/eps;

blue   = [0 0.4470 0.7410 0.4];
yellow = [0.9290 0.6940 0.1250 0.4];
purple = [0.4940 0.1840 0.5560];
green = [96 181 53]/255;

open('Results/approxMFsols.fig');
hold on;
plotcriticalmanifold;
hold off;
xlabel('$K$','Interpreter','LaTeX');
ylabel('$V$','Interpreter','LaTeX');
zlabel('$S$','Interpreter','LaTeX');
interval=[-0.3 0.05 -3 1 0 0.6];
view([50 22]);
axis(interval);

hold on;
idDown = 27; 
fileName = sprintf('RasterData/rasterData-%06i.mat',idDown);
load(fileName);
fprintf('%.22f\n',A0)
IApp = @(t) A0*sin(eps*t);
K = (etaBar + IApp(tHist))/sqrt(M); 
plot3(K,vMeanHist,sMeanHist,'LineWidth',1,'color',green); view([50 22]);

idUp = 28; 
fileName = sprintf('RasterData/rasterData-%06i.mat',idUp);
load(fileName);
fprintf('%.22f\n',A0)
IApp = @(t) A0*sin(eps*t);
K = (etaBar + IApp(tHist))/sqrt(M); 
plot3(K,vMeanHist,sMeanHist,'LineWidth',1,'color',purple); view([50 22]);

xticklabels([]); yticklabels([]); zticklabels([]);

% -------------------------------------------------------------------------------------------

fig = figure;
h = plot3(K,vMeanHist,sMeanHist,'color',purple);
% h = scatter3(K,vMeanHist,sMeanHist,0.8,'filled');
%h.MarkerEdgeColor = purple;
%h.MarkerFaceColor = purple;
%h.MarkerFaceAlpha = 0.2;
view([50 22]);
xlabel('$K$','Interpreter','LaTeX');
ylabel('$V$','Interpreter','LaTeX');
zlabel('$S$','Interpreter','LaTeX');
hold on;
plotcriticalmanifold;
% view([50 22]);
% axis(interval);
interval=[-1 0 -6 6 0 8];
axis(interval);
grid on;
