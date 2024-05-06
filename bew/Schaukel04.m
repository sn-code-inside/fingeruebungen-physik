% -------------------------------------------------------------------------
% Schaukel04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Parametrischer Oszillator
% 
% Lösung der DGL für den Parametrischen Oszillator
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Erwachsener
omega0  = 1;                        % Eigenfrequenz 
kappa   = 0.2;                      % Kopplung
domega  = 0.01;                     % Mismatch
tmax    = 100;
NPts    = 500;
t       = linspace(0,tmax,NPts);    % Zeitbereich
step    = tmax/NPts;
TSPAN   = 0:step:tmax;
kend    = 3;
A0      = 0.01;


% Analytische Näherung
figure()
hold on
for k=1:kend
  domega  = 0.01+(k-1)*0.02;        % Mismatch
  Omega   = omega0+domega/2;
  lambda  = 0.5*sqrt(omega0^2*kappa^2/4 - domega^2);
  check   = omega0*kappa/2/domega;
  x  = A0*exp(lambda.*t).*cos(Omega.*t);
  h(k)= plot(t,x,'Color',Colors(k+1,:),'LineWidth',1);
  strtemp      = num2str(domega/omega0,4);   
  lgdstr(k,:)  = string(sprintf('\\Delta\\omega/\\omega_0 %s ', strtemp));
  if k == 1 
      env = A0*exp(lambda.*t);
      plot (t,env, 'Color',Colors(k+1,:),'LineWidth',1);
      plot (t,-env, 'Color',Colors(k+1,:),'LineWidth',1);
  end
end
grid on
legend(h,lgdstr,'location','southwest');
legend box off
xlabel('$t$','interpreter','latex','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',16)
ht=title('Parametrischer Oszillator Näherung');
set(ht,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Numerische Berechnung ungetrieben, ungedämpft
figure(); 
hold on
for k = 1:kend
    domega  = 0.01+(k-1)*0.02;    % Mismatch
    x0      = A0;
    dx0     = 0;  
    % Anfangswerte
    AB=[x0;dx0];       % AB für ode45
    % Numerische Lösung
    options = odeset('AbsTol',1.e-9,'RelTol',1.e-7);
    [T,Y]=ode45(@(t,Y)dgl_para1(t, Y, omega0, kappa, domega), ...
                                               TSPAN, AB,options); 
    h(k)=plot(T,Y(:,1),'Color', Colors(k+1,:), 'LineWidth',1);
end
plot (t,env, 'Color',Colors(2,:),'LineWidth',1);
plot (t,-env, 'Color',Colors(2,:),'LineWidth',1);
grid on
legend(h,lgdstr,'location','southwest');
legend box off
xlabel('$t$','interpreter','latex','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',16)
h=title('Param. Oszillator ungetrieben, ungedämpft');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Numerische Berechnung getrieben, gedämpft
figure();
hold on
for k = 1:3
    domega  = 0.01+(k-1)*0.02;    % Mismatch
    x0      = A0;
    dx0     = 0;  
    kappa   = 1*kappa;
    p1      = omega0/5;
    p3      = 0.1;
    % Anfangswerte
    AB=[x0;dx0];       % AB für ode45
    % Numerische Lösung
    options = odeset('AbsTol',1.e-9,'RelTol',1.e-7);
    [T,Y]=ode45(@(t,Y)dgl_para2(t, Y, omega0, kappa, domega, p1, p3), ...
                                               TSPAN, AB,options); 
    h(k)=plot(T,Y(:,1),'Color', Colors(k+1,:), 'LineWidth',1);
end
plot (t,env, 'Color',Colors(2,:),'LineWidth',1);
plot (t,-env, 'Color',Colors(2,:),'LineWidth',1);

grid on
legend(h,lgdstr,'location','southwest');
legend box off
xlabel('$t$','interpreter','latex','FontSize',16)
ylabel('$x$','interpreter','latex','FontSize',16)
h=title('Param. Oszillator getrieben, Dämpfung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%%
function dY = dgl_para1(t, Y, omega0, kappa, domega)
    % Y(1)-Auslenkung
    % Y(2)-Geschwindigkeit
    dY = [Y(2);...
         -omega0^2*(1+kappa*cos((2*omega0+domega)*t))*Y(1)];
end

function dY = dgl_para2(t, Y, omega0, kappa, domega, p1, p3)
    % Y(1)-Auslenkung
    % Y(2)-Geschwindigkeit
    dY = [Y(2);...
         -p1*Y(2)-omega0^2*(1+kappa*cos((2*omega0+domega)*t))*Y(1)+ ...
          p3*cos((1*omega0+1*domega)*t)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

