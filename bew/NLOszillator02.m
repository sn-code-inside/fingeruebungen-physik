% -------------------------------------------------------------------------
% NLOszillator02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Getriebener Nichtlinearer Oszillator/Pendel
% 
% Ljapunow-Exponenten beim NLO
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

NPer  = 50;
TA    = 1;
omegA = 2*pi/TA;
omeg0 = 1.25*omegA;
gamma = omeg0/5;
Gamma = 1;              % eine Trajektorie
% Gamma = 1.27270;      % zwei Trajektorien
% Gamma = 1.28;         % zwei Trajektorien
% Gamma = 1.2941525;    % vier Trajektorien
% Gamma = 1.295;        % viele Trajektorien
% Gamma = 1.3;          % wieder stabiler
% Gamma = 1.5;          % wieder stabil
% Gamma = 1.7;          % Chaos
x0    = pi/3;
x0dot = 0;

NPts = 5000;
tmax  = NPer*TA;
NStart = 1000;
tu = 10;
delta0 = 5e-6;
AB1=[x0;x0dot];             % AB für ode45
AB2=[x0+delta0;x0dot];      % AB für ode45
tspan = 0.0:tmax/NPts:tmax;
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);      
[t,x1]=ode45(@dgl_pendel,tspan, AB1,opt,omeg0,gamma,Gamma,omegA); %Numerische Lösung
[t,x2]=ode45(@dgl_pendel,tspan ,AB2,opt,omeg0,gamma,Gamma,omegA); %Numerische Lösung
minx=min(min(x1(:,1),x2(:,1)));
maxx=max(max(x1(:,1),x2(:,1)));
miny=min(min(x1(:,2),x2(:,2)));
maxy=max(max(x1(:,2),x2(:,2)));

figure();
subplot(3,1,1)
% Trajektorien
%hp(1) =  plot(t,wrapToPi(x1(:,1)),'Color',Colors(2,:), 'LineWidth',1);
% hold on
%hp(1) =  plot(t,wrapToPi(x2(:,1)),'Color',Colors(4,:), 'LineWidth',1);
% axis([0 tmax -pi +pi])
hp(1) = plot(t,x1(:,1),'Color',Colors(2,:), 'LineWidth',1);
hold on
hp(2) = plot(t,x2(:,1),'Color',Colors(4,:), 'LineWidth',1);
grid on
axis([tu tmax minx maxx])
ylabel('$q$ in rad','interpreter','latex','FontSize',14)
xlabel('$t$ in s','interpreter','latex','FontSize',14)
hl = legend(hp,' $q_1$' ,' $q_2$','interpreter','latex'); 
set(hl,'FontSize',12);
set(hl,'location', 'northeast');
legend box off;
Gammastr = cat(2,'\Gamma = ', num2str(Gamma,5));
ht = text(tu*1.1, minx + (maxx-minx)/8, Gammastr); 
set(ht,'FontSize',12);
grid on


% Ljapunov Exponent
logDeltat = log(abs(x1(:,1) - x2(:,1)));
lambda = (logDeltat - log(delta0))./t;
y = (logDeltat(NStart:NPts) - log(delta0)); 
x = t(NStart:NPts); %Accidents per state
X = [ones(length(x),1) x];
slopeL =X\y;

% Ljapunov Exponent Nr1
subplot(3,1,2)
plot(t,lambda,'Color',Colors(3,:),'LineWidth',1);
hold on;
miny  = min(min(lambda(NStart:NPts)),-0.1);
maxy  = max(max(lambda(NStart:NPts)),0.1);
axis([tu tmax miny maxy])
line([tu tmax],[0 0],'color',Colors(2,:),'LineStyle',Style(1),...
           'Linewidth',1);
hla1 = line([tu tmax],[lambda(NPts) lambda(NPts)],'color',Colors(4,:),...
           'LineStyle', Style(3),'Linewidth',2);
xlabel('$t$ in s','interpreter','latex','FontSize',14)
ylabel('$\lambda $','interpreter','latex','FontSize',14)
grid on;
lambdastr = cat(2,'$\lambda_1$ = ', num2str(lambda(NPts),5));
hw1 = legend(hla1, lambdastr ,'interpreter','latex'); 
set(hw1,'FontSize',12);
set(hw1,'location', 'best');
legend box off;

% Ljapunov Exponent Nr2
subplot(3,1,3)
plot(t,logDeltat,'Color',Colors(3,:),'LineWidth',1);
hold on;
miny  = min(min(logDeltat(NStart:NPts)),-5);
maxy  = max(max(logDeltat(NStart:NPts)),5);

line([tu tmax],[0 0],'color',Colors(2,:),'LineStyle',Style(1),...
           'Linewidth',1);
hla2= plot(x,+1*log(delta0)+slopeL(1)+slopeL(2)*x,'color',Colors(4,:),...
           'LineStyle',Style(3),'Linewidth',2);
xlabel('$t$ in s','interpreter','latex','FontSize',14)
ylabel('$\ln{|\Delta_{qt}|}$','interpreter','latex','FontSize',14)
axis([tu tmax miny maxy])
grid on;
lambdastr = cat(2,'$\lambda_2$ = ', num2str(slopeL(2),5));
hw2 = legend(hla2, lambdastr ,'interpreter','latex'); 
set(hw2,'FontSize',12);
set(hw2,'location', 'best');
legend box off;
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
