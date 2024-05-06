% -------------------------------------------------------------------------
% Rollenpendel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung des Rollenpendels für
% verschiedene Parameter (ohne Reibung). 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;


% Anfangsbedingungen, Parameter
L=2;                                % Länge Pendel
m1=1; m2=1; M=m1+m2; 
gamma= m2/M;                        % Massen
g=9.81;                             % Schwerebeschleunigung

omega0 = sqrt(g/L);                 % Kreisfrequenz Normalpendel
TP = 2*pi*sqrt(L/g);                % Schwingungsperiode Normalpendel
tmax=3*TP; tstep= TP/50;            % Gesamtzeit und Schrittweite

tsim =(0.0:tstep:tmax);
N=length(tsim);                     % Berechnungszeitraum

x10=0; x1_dot0=-0.2;                % AB Koordinate x1

phi0=-60;                           % Anfangswinkel
if phi0 > 90 
    phi0 = 90; 
end
phi0=deg2rad(phi0); 
phi_dot0=0.25;                      % Anfangsgeschwindigkeit in rad/s

AB =[phi0;phi_dot0];                % Anfangsbedingungen
AB0=[deg2rad(10);phi_dot0];         % AB Normalpendel (kleiner Winkel)

% MATLAB Runge-Kutta (4,5) Ode45 Verfahren

% Default Tolerance     
[tout1,yout1]=ode45(@Diff_Equ,[0 tmax],AB,[],omega0,g,L,gamma);  
% Default Tolerance und ruhendes Pendel    
gamma = 0;
[tout2,yout2]=ode45(@Diff_Equ,[0 tmax],AB0,[],omega0,g,L,gamma);  

% Ausgabe des Pendelwinkels über die Zeit
figure()
str1=cat(2,'Rollenpendel: \phi_0 = ',num2str(rad2deg(phi0),3),...
        '°,  (d\phi/dt)_0 = ',num2str(phi_dot0,3),...
        ' rad/s,  x_1(0) = ', num2str(x10,2),... 
        ' m,  (dx_1/dt)_0) = ',num2str(x1_dot0,3.1),' m/s');
str2=cat(2,'L = ',num2str(L,3),' m , ');
str3=cat(2,'  m_1 = ',num2str(m1,3), ' kg ,  m_2 = ',num2str(m2,3),' kg ',...
           ',  g = ',num2str(g,3),' m/s^2'); 
subplot(3,1,1), 
plot(tout1,rad2deg(yout1(:,1)),'Color', Colors(2,:),'LineStyle','-', ...
    'LineWidth',2);
hold on
plot(tout2,rad2deg(yout2(:,1)*phi0/deg2rad(10)),'Color', Colors(3,:),...
    'LineStyle',':','LineWidth',1); 
xlabel ('t in s','FontSize',12), ylabel('\phi in °','FontSize',12)
title(str1,'FontSize',12)
ym=rad2deg(min([yout1(:,1);yout2(:,1)]));
yp=rad2deg(max([yout1(:,1);yout1(:,1)]));    
axis([0,tmax,ym*(1+0.4),yp*(1+0.4)])
text(.2,ym*(1+0.15),strcat(str2,str3));
h=legend('Rollenpendel','Normalpendel','location','southeast'); 
set(h,'FontSize',12)
legend box off;

% Erhaltungsgröße Px
Px = M*x1_dot0 + m2*L*phi_dot0*cos(phi0);

% Koordinate x1
x1 = x10 + Px*tout1/M - (m2*L/M)*(sin(yout1(:,1))-sin(phi0));
subplot(3,1,2), plot(tout1,x1,'Color', Colors(2,:),'LineStyle','-',...
    'LineWidth',1);
xlabel('t in s','FontSize',12);
ylabel('x_1 in m','FontSize',12);
xlim([0,tmax]);

% Position Pendel x2, z2
x2 = x1  + L*sin(yout1(:,1));
z2 = -L*cos(yout1(:,1));
subplot(3,1,3), plot(x2,z2,'Color', Colors(2,:),'LineStyle','-',...
    'LineWidth',1);
xlabel('x_2 in m','FontSize',12);
ylabel('z_2 in m','FontSize',12);
ylim([-1.1*L,0]);

 
%%
%--------------------------------------------------------------------------
% Simulation 

% Koordinaten über Zeit
z1 = 0*tout1;

% Achsenbestimmung
am=min([x1;x2]); % window size
ap=max([x1;x2]); % window size
if am > -1.1*L amp = -1.1*L; end
if ap  -1.1*L amp = -1.1*L; end
am = floor(am);
ap = round(ap);
figure();
axis([am,ap,-1.2*L,0])
hold on
grid on
pbaspect([(ap-am)/1.2/L 1 1])                   % Aspektverhältnis Achsen
% Aufhängepunkt
plot(x10,0,'o','MarkerSize',10, 'color', Colors(3,:),'LineWidth',2);                       
% Pendelarm 
line([x1(1),x1(1)],[0,-L],'color', Colors(2,:),'LineWidth',2);
% Masse m2 at x2,z2
p(1)=plot(x1(1),-L,'d','color', Colors(2,:),'LineWidth',2);                       
title(str1,'FontSize',12), 
xlabel('x_1, x_2 in m','FontSize',12), ylabel('z in m','FontSize',12)
pause(1.5)

for k=1:length(tout1) 
   clf
   axis([am,ap,-1.2*L,0])
  hold on
   grid on
   pbaspect([(ap-am)/1.2/L 1 1])                    %A spektverhältnis Achsen
   % Aufhängepunkt
   plot(x1(k),0,'o','MarkerSize',10, 'color', Colors(3,:),'LineWidth',2);                       
   % Pendelarm
   line([x1(k),x2(k)],[z1(k),z2(k)],'color', Colors(2,:),'LineWidth',2);
   % Masse m2 at x2,z2
   p(1)=plot(x2(k),z2(k),'d','color', Colors(2,:),'LineWidth',2);                       
   pause(.05)
end
p(2)=plot(x2,z2,'color', Colors(2,:),'LineStyle','-.');   % m2 Trajektorie
pl=legend(p,'m_2', 'Trajektorie','location','northwest');
set(pl,'FontSize',12)
title(str1,'FontSize',12), 
xlabel('x_1, x_2 in m','FontSize',12), ylabel('z in m','FontSize',12)
legend box off
text(am+0.2,-L*1.1,strcat(str2,str3));
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
%--------------------------------------------------------------------------
% Funktionen

% Differentialgleichungssystem
function  dy  =  Diff_Equ(t, y, omega, g, L, gamma)
    dy  =  zeros(2,1);  % es muss ein Spaltenvektor zurückgegeben werden 
    dy(1)  =  y(2);
    dy(2)  = (g/L/(gamma-1))*(1+(gamma*L/g)*y(2)*y(2)*cos(y(1)))*sin(y(1)); 
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

