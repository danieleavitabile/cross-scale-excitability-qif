clear all, close all, clc;

n = 1e4;
list = dir('SimulationData/*.mat');
nFiles = length(list);

A   = zeros(nFiles,1);
solMeas = zeros(nFiles,1);

for i=1:nFiles
  fileName = [list(i).folder '/' list(i).name];
  data = load(fileName);
  A(i) = data.A0;

  % h = histogram(data.spikesHist(:,1),100,'BinLimits',[15 30]);
  % rVals = h.Values/(n*h.BinWidth);
  % solMeas(i) = max(rVals);

  solMeas(i) = max(data.vMeanHist);
end

figure;
plot(A,solMeas,'*-');
xlabel('$A\sqrt{M}$','Interpreter','LaTeX');
ylabel('$\max v$','Interpreter','LaTeX');
