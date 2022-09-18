% Clean
clear all, close all, clc;

%% First dataset with D->D and D->U canard, high K
n      = 1e4;           M = 1e3;
J      = 1/sqrt(M); tauM  = 1;
vt     = 100;        vr   = -vt; 
tauD   = 1.5e-2;   Delta0 = 0.3; 
A0     = 0.51698440*sqrt(M);  eps = 5e-2; %0.5;
etaBar   = -0.5*sqrt(M);
DeltaEta = 1e-4*sqrt(M);
idPerm   = randperm(n)';

% D-->D  0.5; 0.5005; 0.5001; 0.50125; 0.50135; 0.50137; 0.501375; 0.50138;
%        0.5014; 0.50145; 0.50146; 0.501465; 0.501466; 0.5014661; 0.5014662;
%        0.5014663; 0.5014664; 0.501466425; 0.5014664275; 0.501466429; 
%        0.501466432; 0.501466434; 0.50146643475; 0.50146643495; 0.501466434975;
%        0.501466434995; 0.5014664349975; 0.5014664349995; 0.50146643499975;
%        0.501466435;0.501466440; 0.501466442; 0.5014664425;
%        
% D-->U 0.55; 0.51; 0.505; 0.502; 0.5015; 0.50148; 0.5014675; 0.50146725;
%       0.5014672; 0.5014671; 0.50146705; 0.501469; 0.5014675; 0.501467;
%       0.5014665; 0.501466445; 0.501466443; 0.50146644275;

% A0Vals = ...
% [...
%  0.50146644275;...
%         0.5014663; 0.501466432; 0.501466434995; 0.501466435;0.501466440; 0.501466442;
% 	0.5014664425;...
% ]*sqrt(M);

A0Vals = 0.50146644275;
A0Vals = A0Vals*sqrt(M);

rasterFig = figure(); rasterAx = [];
voltFig   = figure();
rateFig   = figure();
sFig      = figure();
pplaneFig = figure();
kvsFig = figure();
kvrFig    = figure();

blue   = [0 0.4470 0.7410 0.4];
yellow = [0.9290 0.6940 0.1250 0.4];
purple = [0.4940 0.1840 0.5560 0.4];

nA0Vals = length(A0Vals);

for ii = 1:nA0Vals

  A0 = A0Vals(ii);

  % Applied current (and derivative)
  IApp = @(t) A0*sin(eps*t); IAppDot = @(t) eps*A0*cos(eps*t);

  % Initial conditions
  t = 0; v = -5.6*ones(n,1); s = 0.43*ones(n,1); 

  % Numerical parameters
  T = 2*pi/eps;
  dt = 1e-3; numSteps = ceil(T/dt); numOut = 1; numPrint = 10000; 
  saveSpikes = false; generateRandPars = false;

  % Random parameters
  if generateRandPars
    E   = ConnectivityMatrix(n,Delta0,M);
    eta = etaBar + DeltaEta*tan(pi*(rand(n,1) - 0.5));
    save('randomPars.mat','E','eta');
    return;
  else
    pars = load('randomPars_n_1e4_M_1e3_etabar_minus0.5.mat');
    % pars = load('randomPars.mat');
    E   = pars.E;
    eta = pars.eta;
  end
   
  % History
  numTimeOut = 1 + floor(numSteps/numOut); 
  itOut = 1;
  tHist     = zeros(numTimeOut,1); tHist(itOut)     = t;
  vMeanHist = zeros(numTimeOut,1); vMeanHist(itOut) = mean(v);
  sMeanHist = zeros(numTimeOut,1); sMeanHist(itOut) = mean(s)/M;
  rateHist  = zeros(numTimeOut,1); rateHist(itOut)  = 0;
  if saveSpikes
    spikesHist = [];
  end

  % Main for loop
  for it = 1:numSteps

    % Euler step
    v = v + dt * ( (eta + v.^2 + IApp(t))/tauM + J*s);
    s = s*(1- dt/tauD);

    % Computing mean
    vMean = mean(v);
    sMean = mean(s)/M;

    % Find id of neurons that fired
    idSpikes = find( v > vt );
    numSpikes  = length(idSpikes);
    
    % Reset
    if (numSpikes > 0)
      v(idSpikes) = vr;
      s = s + sum( E(:,idSpikes) ,2)/tauD;
    end

    % Book keeping
    t = t + dt;

    % Print on screen
    if (mod(it,numPrint) == 0)
      fprintf('%0.3f %0.3f %0.3f\n',t,vMean,sMean);
    end

    % Record in history
    if (mod(it,numOut) == 0)
      itOut            = itOut + 1;
      tHist(itOut)     = t;
      vMeanHist(itOut) = vMean;
      sMeanHist(itOut) = sMean;
      rateHist(itOut)  = numSpikes/(n*dt);
      if saveSpikes
	spikesHist = [spikesHist; t*ones(numSpikes,1) idPerm(idSpikes)];
      end
    end

  end

  % Plot raster plot
  if (saveSpikes)
    figure(rasterFig);
    ax = subplot(nA0Vals,1,ii);
    rasterAx = [rasterAx ax];
    scatter(spikesHist(:,1),spikesHist(:,2),'.'); 
    axis tight;
    xlim([0 t]); ylim([0 n]);
    xlabel('$t$','Interpreter','LaTeX');
    ylabel('$n$','Interpreter','LaTeX','rot',0);
  end

  % Plot means
  figure(voltFig); hold on;
  % plot(tHist,[vMeanHist, sMeanHist, IApp(tHist), rateHist]);
  % plot(tHist,[vMeanHist, sMeanHist]);
  plot(tHist,vMeanHist,'color',blue); hold on; plot(tHist,0*tHist,'k'); hold off;
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  legend({'vMean'},'Interpreter','LaTeX');
  hold off;

  figure(rateFig); hold on;
  % plot(tHist,[vMeanHist, sMeanHist, IApp(tHist), rateHist]);
  % plot(tHist,[vMeanHist, sMeanHist]);
  plot(tHist,rateHist,'.','color',yellow);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$R$','Interpreter','LaTeX');
  legend({'r(t)'},'Interpreter','LaTeX');
  hold off;

  figure(sFig); hold on;
  % plot(tHist,[vMeanHist, sMeanHist, IApp(tHist), rateHist]);
  % plot(tHist,[vMeanHist, sMeanHist]);
  plot(tHist,sMeanHist,'.','color',purple);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$S$','Interpreter','LaTeX');
  legend({'s(t)'},'Interpreter','LaTeX');
  hold off;

  figure(pplaneFig); hold on;
  idPlot    = 1:length(tHist); %12000:length(tHist); 
  tHist     = tHist(idPlot); 
%   vMeanHist = vMeanHist(idPlot); 
  sMeanHist = sMeanHist(idPlot); 
  % rateHist  = rateHist(idPlot);
  % 
  K = (etaBar + IApp(tHist))/sqrt(M); 
  q = IAppDot(tHist)/eps;
  plot3(K,q,sMeanHist,'color',blue); %view([-1 0]);
  hold off;


  figure(kvsFig);hold on;
  plot3(K,vMeanHist,sMeanHist,'color',blue); view([50 22]);
  hold off;
  xlabel('$K$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  zlabel('$S$','Interpreter','LaTeX');

%   figure(kvrFig); hold on;
%   plot3(K,vMeanHist,rateHist,'color',blue); %view([-1 0]);
%   hold off;
%   xlabel('$K$','Interpreter','LaTeX');
%   ylabel('$V$','Interpreter','LaTeX');
%   zlabel('$R$','Interpreter','LaTeX');

  % pause

end
figure(kvsFig), hold on;
plotcriticalmanifold;
interval=[-1 1 -6 8 0 1];
axis(interval);

% figure(rasterFig); linkaxes(rasterAx,'xy');


