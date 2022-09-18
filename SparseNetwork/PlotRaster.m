function PlotRaster(id,flag)

  fileName = sprintf('RasterData/rasterData-%06i.mat',id);
  load(fileName);

  purple = [0.4940 0.1840 0.5560];
  green = [96 181 53]/255;

  switch flag
    case 'up'
      col = purple;
      binLim = [15 45];
      numBins = 200;
      yLim = [0 8000];
    case 'down'
      col = green;
      binLim = [15 30];
      numBins = 100;
      yLim = [0 500];
    otherwise
      error('unknown flag');
  end

  eps = 0.1;
  n      = 1e4; 
  T = 2*pi/eps;


  rasterFig = figure();
  % ax = subplot(nA0Vals,1,ii);
  % rasterAx = [rasterAx ax];
  scatter(spikesHist(:,1),spikesHist(:,2),0.05,'.','MarkerEdgeColor',col,'MarkerFaceColor',col); 
  axis tight;
  xlim([0 T]); ylim([0 n]);
  % xlabel('$t$','Interpreter','LaTeX');
  % ylabel('$n$','Interpreter','LaTeX','rot',0);
  %title(sprintf('$A=%15e$',A0),'Interpreter','LaTeX');
  yticks([]);
  xticks([]);
  box on;

  histFig = figure();
  h=histogram(spikesHist(:,1),numBins,'BinLimits',binLim,'FaceColor',col,'EdgeColor','none');
  % xlabel('$t$','Interpreter','LaTeX');
  %title(sprintf('$A=%15e$',A0),'Interpreter','LaTeX');
  xlim([0 T]);
  ylim(yLim);
  xticks([]);
  yticks([]);
  box on;
  max(h.Values)
end
