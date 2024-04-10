% -------------------------------------------------------------------------
% InvertiertesPendel02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Invertiertes Pendel
% ohne Reibung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
g     = 9.81;                            % g
L     = 1;                               % Länge
q0    = 0.01;                            % Auslemkung t = 0 
qdot0 = 0.0;                             % Anfangsgeschwindigkeit
omega0= sqrt(g/L);                       % Eigenfrequenz
AB    =[q0,qdot0];                       % Anfangsbedingungen
A     = 0.005; 
a     = A/L;                             % Amplitude extern
omegaA= 100;

E0    = qdot0^2+omega0^2*q0^2;
tmax  = 15*2*pi/omega0;                   % Maximalzeit
NPkt  = 1000;                             % Anzahl Punkte
step  = tmax/NPkt;
TSPAN = 0:step:tmax;
kend  = 1;                                               
omegaG = sqrt(2*omega0^2/a^2);            % Grenzfrequenz

%% 
% MATLAB Runge-Kutta (4,5) Ode45 Verfahren
% Fall ohne Dämpfung
    
for k=1:kend
  omega(k) = 10*k*omegaA;                      
  [tout,y] = ode45(@Diff_Equ0,TSPAN,AB,[],omega0, a, omega(k));  
  t0(:,k)     = tout(:);
  theta0(:,k) = y(:,1);   
%   %Kleinwinkelnäherung
%   [tout,y] = ode45(@Diff_Equ0KW,TSPAN,AB,[],omega0, a, omega(k));  
%   t0(:,k)     = tout(:);
%   theta0KW(:,k) = y(:,1);   
end
Omega_eff = sqrt(a^2*omega(kend)^2/2-g/L);
C0 = q0/(1+a);
thetaA=C0*cos(Omega_eff*TSPAN).*(1+a*cos(omega(kend)*TSPAN));

%MATLAB Runge-Kutta (4,5) Ode45 Verfahren
%mit Reibung (Stokes)  
tmaxS  =2*tmax;                         %längere Maximalzeit
NPkt  = 2*1000;                         %Anzahl Punkte
step  = tmaxS/NPkt;
TSPAN2 = 0:step:tmaxS;
for k=1:kend
  omega(k) = 10*k*omegaA;                      
  eta = omega0/8;
  [tout,y] = ode45(@Diff_EquS,TSPAN2,AB,[],omega0, a, omega(k), eta);  
  tS(:,k)     = tout(:);
  thetaS(:,k) = y(:,1);   
end
max_plot = 0.002;

% Aufbereiten für Graphik
for k=1:kend
  for m=1:length(thetaS(:,k))
  if abs(thetaS(end-m,k)) > deg2rad(max_plot) 
      mbegin = m;
      break;  
  end
  end
  for m=1:mbegin
    thetaS(m,k) = NaN;
  end
end
               
fprintf('\n ');
fprintf('\n L   = %8.3f m', L);
fprintf('\n a   = %8.3f ', a);
fprintf('\n ');
fprintf('\n omega0^2  = %8.2f 1/s^2', omega0^2);
fprintf('\n eta       = %8.2f 1/s', eta);
fprintf('\n omega     = %8.2f 1/s', omega(kend));
fprintf('\n Omega_eff = %8.2f 1/s', Omega_eff);
fprintf('\n T_eff     = %8.2f s', 2*pi/Omega_eff);
fprintf('\n ');
        
                                                  
%%
% Graphische Ausgabe

figure();
% ohen Reibung
subplot(1,2,1)
hold on
for k = 1:kend
 plot(t0(:,k),rad2deg(theta0(:,k)), 'Color',Colors(k+1,:),'Linewidth',1,...
     'LineStyle',Style(1));
 lgdstr(k,:)= strcat('\omega =', ...
             string(num2str(omega(k),'% 8.0f 1/s')));
%  plot(t0(:,k),rad2deg(theta0KW(:,k)), 'Color',Colors(2*k+1,:),'Linewidth',1,...
%      'LineStyle',Style(3));
%  lgdstr(2*k,:)= "Kleinwinkelnäherung";
end
plot(TSPAN,rad2deg(thetaA), 'Color',Colors(4,:),'Linewidth',1,...
     'LineStyle',Style(1));
lgdstr(kend+1,:)= "analyt. Näherung";
ylabel('Auslenkung  \theta in ° bzw. rad','FontSize',14);
xlabel('Zeit t in s ','FontSize',14);
legend(lgdstr,'location','northwest');
legend box off
axis([0 tmax -inf inf])
grid on
title('Inverses Pendel','FontWeight','normal',...
      'FontSize',14);
set(gca,'Fontsize', 14);

% mit Reibung
subplot(1,2,2)
hold on
for k = 1:kend
 plot(tS(:,k),rad2deg(thetaS(:,k)), 'Color',Colors(k+1,:),'Linewidth',1,...
     'LineStyle',Style(1));
 lgdstr(k,:)= strcat('\omega =', ...
             string(num2str(omega(k),'% 8.0f 1/s')));
end
ylabel('Auslenkung  \theta in ° bzw. rad','FontSize',14);
xlabel('Zeit t in s ','FontSize',14);
legend(lgdstr(1:kend),'location','northwest');
legend box off
axis([tmax tmaxS -0.005 0.005])
grid on
title('Inverses Pendel','FontWeight','normal',...
      'FontSize',14);
set(gca,'Fontsize', 14);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%%
% Differentialgleichungssysteme

% exakt ohne Reibung
function  dy  =  Diff_Equ0(t, y, omega0, a, omega)
    dy     =  zeros(2,1); % es muss ein Spaltenvektor zurückgegeben werden 
    dy(1)  =  y(2);
    dy(2)  = -sin(y(1))*(a*omega^2*cos(omega*t)-omega0^2); 
end

% numerisch Kleinwinkelnäherung
function  dy  =  Diff_Equ0KW(t, y, omega0, a, omega)
    dy     =  zeros(2,1); % es muss ein Spaltenvektor zurückgegeben werden
    dy(1)  =  y(2);
    dy(2)  = -y(1)*(a*omega^2*cos(omega*t)-omega0^2); 
end

% analytische Näherung
function  dy  =  Diff_EquS(t, y, omega0, a, omega, eta)
    dy     =  zeros(2,1); % es muss ein Spaltenvektor zurückgegeben werden 
    dy(1)  =  y(2);
    dy(2)  = -sin(y(1))*(a*omega^2*cos(omega*t)-omega0^2) - eta*y(2); 
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


