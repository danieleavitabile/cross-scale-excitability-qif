clear all, close all, clc;
idVals = [4 10 27 28];
for id = idVals
  if id > 27 
    PlotRaster(id,'up');
  else
    PlotRaster(id,'down');
  end
end
