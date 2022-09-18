clear all, close all, clc;

blue   = [0 0.4470 0.7410 0.4];
yellow = [0.9290 0.6940 0.1250 0.4];
purple = [0.4940 0.1840 0.5560 0.4];
green = [96 181 53]/255;

n = 1e4;
list = dir('SimulationData/*.mat');
nFiles = length(list);

eps = 0.1;
n      = 1e4; 
T = 2*pi/eps;

A   = zeros(nFiles,1);
solMeas = zeros(nFiles,1);

figure, hold on;
for i=1:28
  fileName = [list(i).folder '/' list(i).name];
  data = load(fileName);
  A(i) = data.A0;
  fprintf('%.22f\n',A(i))

  % h = histogram(data.spikesHist(:,1),100,'BinLimits',[15 30]);
  % rVals = h.Values/(n*h.BinWidth);
  % solMeas(i) = max(rVals);

  %solMeas(i) = max(data.vMeanHist);
  if i < 28
    plot(data.tHist,data.vMeanHist,'color',green);
  else
    plot(data.tHist,data.vMeanHist,'color',purple);
  end
end
hold off;
xlabel('$t$','Interpreter','LaTeX');
ylabel('$v$','Interpreter','LaTeX');
box on;
xlim([0 T]);
