% Clean
clear all, close all, clc;

%% First dataset with D->D and D->U canard, high K
n      = 1e4;           M = 1e3;
J      = 1/sqrt(M); tauM  = 1;
vt     = 100;        vr   = -vt; 
tauD   = 1.5e-2;   Delta0 = 0.3; 
eps = 0.1; %0.5;
etaBar   = -0.5*sqrt(M);
DeltaEta = 1e-4*sqrt(M);
idPerm   = randperm(n)';

data1 = load('A0Vals-1.mat');
data2 = load('A0Vals-2.mat');
A0Vals = sort([data1.A0Vals; data2.A0Vals]);
% A0Vals = [A0Vals(1); A0Vals(end)];

nA0Vals = length(A0Vals);
intSVals = zeros(size(A0Vals));

% rasterFig = figure(); rasterAx = [];
voltFig   = figure();
rateFig   = figure();
sFig      = figure();
% pplaneFig = figure();
kvsFig = figure();
% kvrFig    = figure();
bd = figure();

blue   = [0 0.4470 0.7410 0.4];
yellow = [0.9290 0.6940 0.1250 0.4];
purple = [0.4940 0.1840 0.5560 0.4];

for ii = 1:nA0Vals

  %format long e
  A0 = A0Vals(ii);

  % Applied current (and derivative)
  IApp = @(t) A0*sin(eps*t); IAppDot = @(t) eps*A0*cos(eps*t);

  % Initial conditions
  t = 0; v = -5.6*ones(n,1); s = 0.43*ones(n,1); 

  % Numerical parameters
  % T = 2*pi/eps;
  T = 2*pi/eps;
  dt = 1e-4; numSteps = ceil(T/dt); numOut = 1; numPrint = 100000; 
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
  plot(tHist,vMeanHist,'color',blue); hold on; plot(tHist,0*tHist,'k'); hold off;
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  legend({'vMean'},'Interpreter','LaTeX');
  hold off;

  figure(rateFig); hold on;
  plot(tHist,rateHist,'.','color',yellow);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$R$','Interpreter','LaTeX');
  legend({'r(t)'},'Interpreter','LaTeX');
  hold off;

  figure(sFig); hold on;
  plot(tHist,sMeanHist,'.','color',purple);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$S$','Interpreter','LaTeX');
  legend({'s(t)'},'Interpreter','LaTeX');
  hold off;

 %  figure(pplaneFig); %hold on;
 %  idPlot    = 1:length(tHist); %12000:length(tHist); 
 %  tHist     = tHist(idPlot); 
%%    vMeanHist = vMeanHist(idPlot); 
 %  sMeanHist = sMeanHist(idPlot); 
 %  % rateHist  = rateHist(idPlot);
 %  % 
 K = (etaBar + IApp(tHist))/sqrt(M); 
 %  q = IAppDot(tHist)/eps;
 %  plot3(K,q,sMeanHist,'color',blue); %view([-1 0]);
 %  hold off;

  figure(kvsFig); hold on;
  plot3(K,vMeanHist,sMeanHist,'color',blue); view([50 22]);
  hold on;
  plotcriticalmanifold;
  hold off;
  interval=[-0.5 0.5 -3 0.5 0 0.6];
  axis(interval);
  xlabel('$K$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  zlabel('$S$','Interpreter','LaTeX');
  axis(interval);


%   figure(kvrFig); hold on;
%   plot3(K,vMeanHist,rateHist,'color',blue); %view([-1 0]);
%   hold off;
%   xlabel('$K$','Interpreter','LaTeX');
%   ylabel('$V$','Interpreter','LaTeX');
%   zlabel('$R$','Interpreter','LaTeX');

  figure(bd);
  intSVals(ii) = sum(sMeanHist)*dt;
  plot(A0Vals(1:ii),intSVals(1:ii)/T,'.-','Color',blue);
  xlabel('$A_0$','Interpreter','LaTeX');
  ylabel('$\Vert s \Vert_1$','Interpreter','LaTeX');

  % pause

end

% figure(rasterFig); linkaxes(rasterAx,'xy');
