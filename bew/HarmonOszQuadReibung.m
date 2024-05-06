% -------------------------------------------------------------------------
% HarmonOszQuadReibung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Beispielberechnungen zum Harmonischen Oszillator
% mit linearer und quadratischer Reibung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%Parameter
m     = 0.2;                             %Masse
k     = 1;                               %rücktreibende Kraftkonstante
q0    = 4.0;                             %Auslemkung t = 0 
qdot0 = 0.0;                             %Anfangsgeschwindigkeit
omega0= sqrt(k/m);                       %Eigenfrequenz
AB    =[q0,qdot0];                       %Anfangsbedingungen
E0    = qdot0^2+omega0^2*q0^2;
tmax  = 8;                               %Maximalzeit
NPkt  = 500;                             %Anzahl Punkte
t     = linspace(0,tmax,NPkt);           %Zeitbereich
step  = tmax/NPkt;
TSPAN = 0:step:tmax;
kend  = 5;                                               
 



%% 
% Fall ohne Dämpfung
eta10 = 0.00;                            %Reibungsterm 1
B     = sqrt(q0^2+qdot0^2/omega0^2);     %Konstante B  
theta = atan(omega0*q0/(qdot0));         %Winkel theta
q     = B.*sin(omega0*t+theta);          %Ungedämpfter Oszillator

                                                 
%% 
% Fall mit Stokes-Dämpfung
    
%MATLAB Runge-Kutta (4,5) Ode45 Verfahren
%Default Tolerance     
for k=1:kend
  etaS(k) = 0.25*(k-1);                         %Reibungsterm 1
  [tout,y] = ode45(@Diff_EquS,TSPAN,AB,[],omega0,etaS(k));  
  tS(:,k)       = tout(:);
  xS(:,k)     = y(:,1);   
  xSdot(:,k)  = y(:,2);
  ES(:,k)     = xSdot(:,k).^2 + omega0^2*xS(:,k).^2;
end

                                                  

%% 
% Fall mit Newton-Dämpfung

AB =[q0,qdot0];                     %Anfangsbedingungen

%MATLAB Runge-Kutta (4,5) Ode45 Verfahren
%Default Tolerance     
for k=1:kend
  etaN(k) = 0.25*(k-1);                         %Reibungsterm 1
  [tout,y] = ode45(@Diff_EquN,TSPAN,AB,[],omega0,etaN(k));  
  tN(:,k)     = tout(:);
  xN(:,k)     = y(:,1);   
  xNdot(:,k)  = y(:,2);
  EN(:,k)     = xNdot(:,k).^2 + omega0^2*xN(:,k).^2;
end
                                                  
%%
% Graphik Auslenkungen

figure();
subplot(1,2,1)
hold on
lgdstr(1,:)= strcat(' \eta_S = ',(num2str(etaS(1),'% 4.2f')));
plot(t,q,'Color',Colors(8,:),'Linewidth',1,'LineStyle',Style(3));
for k = 2:kend
 plot(tS(:,k),xS(:,k), 'Color',Colors(k,:),'Linewidth',1,...
     'LineStyle',Style(1));
 lgdstr(k,:)= strcat(' \eta_S = ',(num2str(etaS(k),'% 4.2f')));
end
ylabel('Auslenkung  \it{q} \rm in m','FontSize',14);
xlabel('Zeit \it{t} \rm in s ','FontSize',14);
legend(lgdstr);
legend box off
axis([0 8 -4 4 ])
grid on
title('Harmonischer Oszillator Stokes-Reibung','FontWeight','normal',...
      'FontSize',14);
set(gca,'Fontsize', 14);

subplot(1,2,2)
lgdstr(1,:)= strcat(' \eta_N = ',(num2str(etaN(k),'% 4.2f')));
plot(t,q,'Color',Colors(8,:),'Linewidth',1,'LineStyle',Style(3));
for k = 2:kend
 hold on
 plot(tN(:,k),xN(:,k), 'Color',Colors(k,:),'Linewidth',1,...
     'LineStyle',Style(1));
 lgdstr(k,:)= strcat(' \eta_N = ',(num2str(etaS(k),'% 4.2f')));
end
ylabel('Auslenkung  \it{q} \rm in m','FontSize',14);
xlabel('Zeit \it{t} \rm in s ','FontSize',14);
legend(lgdstr);
legend box off
axis([0 8 -4 4 ])
title('Harmonischer Oszillator Stokes-Reibung','FontWeight','normal',...
      'FontSize',14);
grid on
set(gca,'Fontsize', 14);

%%
% Graphik Energien

figure();
subplot(1,2,1)
lgdstr(1,:)= strcat(' \eta_S = ',(num2str(etaS(1),'% 4.2f')));
semilogy(tS(:,1),ES(:,1),'Color',Colors(8,:),'Linewidth',1,...
        'LineStyle',Style(3));
for k = 2:kend
 hold on
 semilogy(tS(:,k),ES(:,k), 'Color',Colors(k,:),'Linewidth',1,...
     'LineStyle',Style(1));
 lgdstr(k,:)= strcat(' \eta_S = ',(num2str(etaS(k),'% 4.2f')));
end
ylabel('Energie \it{E} \rm normalisiert','FontSize',14);
xlabel('Zeit \it{t} \rm in s ','FontSize',14);
legend(lgdstr,'location','southwest');
legend box off
axis([0 8 1e-2 100])
grid on
title('Harmonischer Oszillator Stokes-Reibung','FontWeight','normal',...
      'FontSize',14);
set(gca,'Fontsize', 14);

subplot(1,2,2)
lgdstr(1,:)= strcat(' \eta_N = ',(num2str(etaN(1),'% 4.2f')));
semilogy(tS(:,1),EN(:,1),'Color',Colors(8,:),'Linewidth',1,'LineStyle',Style(3));
for k = 2:kend
 hold on
 semilogy(tN(:,k),EN(:,k), 'Color',Colors(k,:),'Linewidth',1,...
     'LineStyle',Style(1));
 lgdstr(k,:)= strcat(' \eta_N = ',(num2str(etaN(k),'% 4.2f')));
end
ylabel('Energie \it{E} \rm normalisiert','FontSize',14);
xlabel('Zeit \it{t} \rm in s ','FontSize',14);
legend(lgdstr,'location','southwest');
legend box off
axis([0 8  1e-2 100])
title('Harmonischer Oszillator Stokes-Reibung','FontWeight','normal',...
      'FontSize',14);
grid on
set(gca,'Fontsize', 14);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% 
% Funktionen
% Differentialgleichungssystem
function  dy  =  Diff_EquN(t, y, omega0, etaN)
    dy  =  zeros(2,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
    dy(1)  =  y(2);
    dy(2)  = -omega0^2*y(1)-etaN*sign(y(2))*y(2)*y(2); 
end

% Differentialgleichungssystem
function  dy  =  Diff_EquS(t, y, omega0, etaS)
    dy  =  zeros(2,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
    dy(1)  =  y(2);
    dy(2)  = -omega0^2*y(1)-etaS*y(2); 
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

