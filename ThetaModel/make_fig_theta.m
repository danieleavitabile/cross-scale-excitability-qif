% Clean
clear all, close all, clc;

J   = 6; tauS = 0.3;
eps = 1e-1; eta = 0.5;

%%%% Eps = 1e-1;

% load('A0Vals.mat');
% A0Vals = [A0Vals; linspace(0.875,0.8875,10)'; linspace(0.8883,0.8887,50)'];
% A0Vals = sort(A0Vals);
% A0Vals = A0Vals(find(0.8875 <= A0Vals & A0Vals <= 0.8906));
% 
% A0Vals
% pause

% A0Vals = linspace(0.8882,0.8892,100); ta = 40; tb = 65;
A0Vals = 0.83; ta = 0; tb = 70;

% A0Vals = [0.1; 3];
% A0Vals = [0.1; 0.6];
% A0Vals = [0.8; 0.9];
intSVals = zeros(size(A0Vals));

blue   = [0 0.4470 0.7410 0.4];
yellow = [0.9290 0.6940 0.1250 0.4];
purple = [0.4940 0.1840 0.5560 0.4];
grey =   [222 223 224]/255;

nA0Vals = length(A0Vals);

vFig   = figure();
sFig   = figure();
insetFig   = figure();
% pplaneFig   = figure();
critmanfig   = figure();
bd = figure();


for ii = 1:nA0Vals

  A0 = A0Vals(ii);

  % Applied current (and derivative)
  IApp = @(t) A0*sin(eps*t); IAppDot = @(t) eps*A0*cos(eps*t);

  % Initial conditions
  t = 0; v = 1.0; s = 0; 

  % Numerical parameters
  % T = 2*pi/eps;
  T = 70;
  tout = t;
  yout = [v s];

  % Plot 
  refine = 4;
  options = odeset('Events',@events,'Refine',refine,'RelTol',1e-9,'AbsTol',1e-9);
  y0 = yout;
  tstart = t;
  tfinal = T;

  while tout(end) < tfinal 

   [t,y,te,~,~] = ode23(@(t,y) f(t,y,IApp,J,tauS,eta),[tstart tfinal],y0,options);

   nt = length(t);
   tout = [tout; t(2:nt)];
   yout = [yout; y(2:nt,:)];

   y0(1) = -pi/2;
   y0(2) = y(nt,2) + 1;

   options = odeset(options,'InitialStep',t(nt)-t(nt-refine),...
      'MaxStep',t(nt)-t(1),'RelTol',1e-9,'AbsTol',1e-9);

    tstart = t(nt);
  end

  intSVals(ii) = yout(1:end-1,2)'*diff(tout);
  id = find(ta <= tout & tout <= tb);

  figure(vFig); hold on;
  plot(tout(id),yout(id,1),'color',blue);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  hold off;

  figure(insetFig); hold on;
  id1 = find(tout > 49.4); id1 = id1(1);
  %id2 = id1+1000;
  dv  = diff(yout(id1:id(end),1)); id2 = find(abs(dv) > 3.1); id2 = id2(1)-1;
  plot(tout(id),yout(id,1),'color',grey);
  plot(tout(id1:id1+id2),yout(id1:id1+id2,1),'color',blue);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  hold off;

  figure(sFig); hold on;
  plot(tout(id),yout(id,2),'color',purple);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$S$','Interpreter','LaTeX');
  hold off;

 %  figure(pplaneFig); hold on;
 %  K = (eta + IApp(tout));
 %  q = IAppDot(tout)/eps;
 %  plot3(K,q,yout(:,2),'color',blue); view([-170 38]);
 %  xlabel('$K$','Interpreter','LaTeX');
 %  ylabel('$Q$','Interpreter','LaTeX');
 %  zlabel('$S$','Interpreter','LaTeX');
 %  hold off;

  figure(critmanfig); hold on;
  id = find(ta <= tout & tout <= tb);
  K = (eta + IApp(tout));
  plot(K(id),yout(id,1),'color',blue);
  if ii == 1
    hold on; 
    vVals = linspace(-1,1,1000);
    KVals = (cos(vVals)-1)./(1+cos(vVals));
    plot(KVals,vVals,'color',[1 0 0]);
    hold off;
  end
  xlabel('$K$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  hold off;

  figure(bd);
  plot(A0Vals(1:ii),1./intSVals(1:ii),'*','Color',blue);
  xlabel('$A_0$','Interpreter','LaTeX');
  ylabel('$\Vert s \Vert_1$','Interpreter','LaTeX');

end


function dydt = f(t,y,IApp,J,tauS,eta)

  theta = y(1);
  s     = y(2);
  dydt = zeros(size(y));
  dydt(1) = 1- cos(theta)+(1+cos(theta))*(eta+IApp(t)+J*s);
  dydt(2) = -s/tauS;

end

function [value,isterminal,direction] = events(t,y)
  value = y(1)-pi/2; 
  isterminal = 1;   
  direction = 1; 
end

