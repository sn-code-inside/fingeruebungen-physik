% -------------------------------------------------------------------------
% NLOszillator01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Getriebener Nichtlinearer Oszillator/Pendel
% 
% Phasenraum/Poincare-Schnitt beim NLO
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

NPer  = 30;
TA    = 1;
omegA = 2*pi/TA;
omeg0 = 1.25*omegA;
gamma = omeg0/5;
Gamma = 1;              % eine Trajektorie
 Gamma = 1.27270;       % zwei Trajektorien
% Gamma = 1.28;         % zwei Trajektorien
 Gamma = 1.2941525;     % vier Trajektorien
% Gamma = 1.295;        % viele Trajektorien
% Gamma = 1.3;          % wieder stabiler
% Gamma = 1.5;          % wieder stabil
% Gamma = 1.7;          % Chaos
x0    = pi/3;
x0dot = 0;

tmax  = NPer*TA;
AB=[x0;x0dot];       % AB für ode45
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-5);      
[t,x]=ode45(@dgl_pendel,[0.0,tmax],AB,opt,omeg0,gamma,Gamma,omegA); %Numerische Lösung
minx=min(x(:,1));
maxx=max(x(:,1));
miny=min(x(:,2));
maxy=max(x(:,2));
figure();
subplot(3,1,1)
plot(t,wrapToPi(x(:,1)),'Color',Colors(2,:), 'LineWidth',2);
hold on
axis([0 tmax -pi pi])
grid on
ylabel('$q$ in rad','interpreter','latex','FontSize',14)
xlabel('$t$ in s','interpreter','latex','FontSize',14)

NLong = length(t)-round(0.25*length(t));
subplot(3,1,2)
plot(x(NLong:end,1),x(NLong:end,2),'Color',Colors(3,:),'LineWidth',2);
xlabel('$q$ in rad','interpreter','latex','FontSize',14)
ylabel('$\dot{q}$ in rad/s','interpreter','latex','FontSize',14)
grid on;

AB1 = [x(end,1),x(end,2)];
NLong = 100*NPer;
tlong = tmax;
tverylong = NLong*TA;
tvector = tlong:(2*pi/omegA):tverylong;
[tPoinc,xPoinc]=ode113(@dgl_pendel,tvector,AB1,opt,omeg0,gamma,Gamma,omegA); %Numerische Lösung
subplot(3,1,3)
plot(wrapToPi(xPoinc(:,1)),xPoinc(:,2),'o','MarkerSize',2,'Color',Colors(4,:),'LineWidth',1);
xlabel('$q$ in rad','interpreter','latex','FontSize',14)
ylabel('$\dot{q}$ in rad/s','interpreter','latex','FontSize',14)
Gammastr = cat(2,'\Gamma = ', num2str(Gamma,5));
if Gamma > 1.5
    axis([-pi pi miny maxy])
    text(-pi+0.2,miny + (maxy-miny)/8, Gammastr);
else
    axis([-pi 0 miny maxy])
    text(-pi+0.2,miny + (maxy-miny)/8, Gammastr);   
end
  
grid on
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
function dY = dgl_pendel(t, Y, omeg0, gamma, Gamma, omegA)
    % Y(1)-Winkel, Y(2)-Winkelgeschwindigkeit
    dY = [Y(2); -omeg0^2*sin(Y(1))-2*gamma*Y(2)+Gamma*omeg0^2*cos(omegA*t)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
