% -------------------------------------------------------------------------
% AsteroidKollision.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Asteroiden-Auftreffwahrscheinlichkeit 
%
% Programm berechnet Trajektorien für Asteroiden mit
% Auftreffwahrscheinlichkeit
% Planet:  Erde, Jupiter
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

% Variablen, Konstanten 
G   = 6.67259e-11;              % Gravitationskonstante in m^3/kg/s^2
mE  = 5.9720e24 ;               % Erdmasse in kg
RE  = 6378000 ;                 % Erdradius in m
GE  = mE*G;                     % spez. Gravitationskonstante in m^3/s^2
mJ  = 1.898e27 ;                % Masse Jupiter in kg
GJ  = mJ*G;                     % spez. Gravitationskonstante in m^3/s^2
RJ  = 71492000 ;                % Jupiterradisu in m
vinf = 14000;                   % initiale Geschwindigkeit m/s

% Erde
kend =5;

theta = linspace(-pi,pi,1000);

figure();
hold on
for k=1:kend
   b   = k*1*RE;                     % Stoßparameter in m
   L   = vinf*b;
   exc = sqrt(1+vinf^4*b^2/GE/GE);
   alpha = acos(-1/exc);
   r   = L^2/GE*(1./(1+exc*cos(theta)))/RE;
   for m=1:length(r) 
       if r(m)  < 0 
           r(m) = NaN;
       end
   end
   x = r.*cos(theta-alpha);
   y = r.*sin(theta-alpha);
   h(k)=plot(x,y, 'color',Colors(k,:),'Linewidth',2);
   strtemp      = num2str(k,2);   
   lgdstr(k,:)  = string(sprintf('b/RE = %s ', strtemp));
end
axis([-5 5 -5 5]);
PlotCircle (0,0,1,Colors(3,:),2);
legend(h,lgdstr,'location','northeast');
legend box off
axis square
grid on
grid minor
% legend(h,lgdstr,'location','northeast');
% legend box off
xlabel('$x$ in $R_E$','interpreter','latex','FontSize',16)
ylabel('$y$ in $R_E$','interpreter','latex','FontSize',16)
vstr = num2str(vinf/1000,3);
vstr = strcat('v_\infty = ',vstr);
vstr1 = strcat(vstr,' km/s');
vstr2 = "Erde  ";
ht=title(strcat(vstr2,vstr1));
set(ht,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

figure();
vinf = 14000;                        % initiale Geschwindigkeit m/s
hold on
for k=1:kend
   b   = k*2*RJ;                     % Stoßparameter in m
   L   = vinf*b;
   exc = sqrt(1+vinf^4*b^2/GJ/GJ);
   alpha = acos(-1/exc);
   r   = L^2/GJ*(1./(1+exc*cos(theta)))/RJ;
   for m=1:length(r) 
       if r(m)  < 0 
           r(m) = NaN;
       end
   end
   x = r.*cos(theta-alpha);
   y = r.*sin(theta-alpha);
   h(k)=plot(x,y,'color',Colors(k,:),'Linewidth',2);
   strtemp      = num2str(k,2);   
   lgdstr(k,:)  = string(sprintf('b/RJ = %s ', strtemp));
end
axis([-5 5 -5 5]);
PlotCircle (0,0,1,Colors(5,:),2)
axis square
grid on
grid minor
% legend(h,lgdstr,'location','northeast');
% legend box off
xlabel('$x$ in $R_J$','interpreter','latex','FontSize',16)
ylabel('$y$ in $R_J$','interpreter','latex','FontSize',16)
vstr = num2str(vinf/1000,3);
vstr = strcat('v_\infty = ',vstr);
vstr1 = strcat(vstr,' km/s');
vstr2 = "Jupiter  ";
ht=title(strcat(vstr2,vstr1));
set(ht,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
