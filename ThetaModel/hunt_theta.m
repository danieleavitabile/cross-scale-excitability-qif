% Clean
clear all, close all, clc;

J   = 6; tauS = 0.3;
eps = 1e-1; eta = 0.5;

% U --> U 
A0ValsUpUp   = 0.8;

% U --> D
A0ValsUpDown = 0.9;

blue   = [0 0.4470 0.7410 0.4];
yellow = [0.9290 0.6940 0.1250 0.4];
purple = [0.4940 0.1840 0.5560 0.4];

pos = [100 900 1000 300];
vFig   = figure('Position',pos);
sFig   = figure('Position',pos);
pplaneFig =  figure('Position',pos);
bd =  figure('Position',pos);

for ii = 1:100

  format long e
  A0 = 0.5*(A0ValsUpDown(end)+A0ValsUpUp(end))

  % Applied current (and derivative)
  IApp = @(t) A0*sin(eps*t); IAppDot = @(t) eps*A0*cos(eps*t);

  % Initial conditions
  t = 0; v = 1.0; s = 0; 

  % Numerical parameters
  % T = 4*pi/eps;
  T = 80;
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

  sVals(ii) = max(yout(:,2))-min(yout(:,2));

  figure(sFig);% hold on;
  plot(tout,yout(:,2),'color',purple);
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  %hold off;

  figure(vFig); %hold on;
  plot(tout,yout(:,1),'color',blue); hold on; plot(tout,0*tout,'r'); hold off;
  xlabel('$t$','Interpreter','LaTeX');
  ylabel('$V$','Interpreter','LaTeX');
  %hold off;

%   figure(pplaneFig); 
%   K = (eta + IApp(tout));
%   q = IAppDot(tout)/eps;
%   % plot(K,yout(:,1),'color',blue); 
%   plot3(K,q,yout(:,1),'color',blue); 
%   xlabel('$K$','Interpreter','LaTeX');
%   ylabel('$Q$','Interpreter','LaTeX');
%   zlabel('$S$','Interpreter','LaTeX');
%  %  hold off;

  figure(bd); hold on;
  plot(A0,1./(yout(1:end-1,2)'*diff(tout)),'*');
  hold off;

 flag = input('Up or Down [u/d]?','s');
 switch flag
   case 'u'
    A0ValsUpUp = [A0ValsUpUp; A0];
   case 'd'
    A0ValsUpDown = [A0ValsUpDown; A0];
   otherwise
     error('unknown flag');
 end

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

