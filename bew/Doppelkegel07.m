% -------------------------------------------------------------------------
% Doppelkegel07.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärtsrollender Doppelkegel
% 
% Programm berechnet inkorrekte nur näherungsweise 
% gültige Lagrange-Funktion und 
% Lagrange-Gleichungen des aufwärtsrollenden Doppelkegels  
% 
% Benutzt die symbolische Euler-Lagrange-Berechnung nach 
% Morten Veng (2021). Euler-Lagrange Solver 
% (https://www.mathworks.com/matlabcentral/fileexchange/93275-euler-lagrange-solver), 
% MATLAB Central File Exchange. Retrieved December 27, 2021.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Definiert die verallgemeinerten Koordinaten theta1, theta2 und Parameter

syms m J R  ya dya A B C D C1 C2 C3 C4 C5 g 'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
Q  = [ya];
dQ = [dya];


% kinetische Energie
T_kin = 1/2*m*dya.*dya.*(A+B+C./(D-ya).^2.*(1 + ya./(D-ya)).^2);

% potentielle Energien
U_pot =  m*g*(C3*R+C5*ya);


% Lagrange-Funktion
L = T_kin - U_pot;

%% Ausgabe der Ergebnisse
fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf("\n");

% Externe Kräfte 
QF = [0];

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(Q,dQ,L,QF,true);


%% Numerik
% Parameter

% Parameter

H       = 6.5;              % Länge Schaukel cm
R       = 3.0;              % Länge Oberkörper cm 
phi     = deg2rad(15.5);    % 
theta   = deg2rad(06.5);    % 
psi     = atan(R/H);        % 
m       = 122.5;
J       = 330.75;
Js      = 3*m*R^2/10;
g       = 981;           % g in cm/s
tmax    = 25;

fprintf('\n ');
fprintf('\n phi   = %4.2f°', rad2deg(phi));
fprintf('\n psi   = %4.2f°', rad2deg(psi));
fprintf('\n theta = % 4.2f°', rad2deg(theta));
fprintf('\n ');
fprintf('\n H  = %8.2f cm', H);
fprintf('\n R  = %8.2f cm', R);
fprintf('\n m  = %8.2f g', m);
fprintf('\n J  = %8.2f gcm²', J);
fprintf('\n g  = %8.2f cm/s', g);
fprintf('\n ');
fprintf('\n Rollbedingung erfüllt ? ');
b1 = tan(theta);
b2 = tan(phi)*tan(psi)/sqrt(1-(tan(phi)*tan(psi))^2);
if b1 < b2
    fprintf('\n %s  < %s ! Rollbedingung erfüllt! ',...
        num2str(b1,4), num2str(b2,4));
    fprintf('\n ');
else
    fprintf('\n %s  > %s ! Rollbedingung nicht erfüllt! ',...
        num2str(b1,4), num2str(b2,4));
    fprintf('\n ');
end

fprintf('\n ');

% Parameter in LGL

C1    = tan(phi)*tan(psi);
C2    = 1+C1*sin(theta);
C3    = cos(theta);
C4    = -R*sin(theta);
C5    = tan(theta)-C1*C3;
A     = C2^2;
B     = C5^2;
C     = J/(m*C3^2*C1^2);
D     = R/C1;

ya0 = 0;
dya0= 0;


%% Berechnungen ODE45

% Anfangswerte
AB=[ya0;dya0]; % AB für ode45


P1.m     = m;
P1.g     = g;
P1.J     = J;
P1.H     = H;
P1.R     = R;
P1.C1    = C1;
P1.C2    = C2;
P1.C3    = C3;
P1.C4    = C4;
P1.C5    = C5;
P1.A     = A;
P1.B     = B;
P1.C     = C;
P1.D     = D;
P1.phi   = phi;

% Vergleich mit numerischer Berechnung aus Doppelkegel02.m (LGL)
QLGL = ImportFile01('Data/DKyzLGL.txt')
tLGLex  =  QLGL.t;
yLGLex  =  QLGL.y;
zLGLex  =  QLGL.z;
dyLGLex =  QLGL.dy;
dzLGLex =  QLGL.dz;


% Numerische Lösung Lagrangegleichung
opt=odeset('AbsTol',1.e-10,'RelTol',1.e-10,'events',@MyEvent);
% Numerische Lösung LGL
[t, Y, TE, YE, IE] = ode45(@dgl_DoppelKegel,[0.0,tmax],AB,opt,P1); 
ya  = Y(:,1);
dya = Y(:,2); 
y = ya*C2+R*sin(theta);
dy = dya*C2;
z = ya*C5+ R*C3;
dz = dya*C5;
kend = length(y);
yr = linspace(0,25);
zr = yr*tan(theta);

ya1  = linspace(0,max(ya),100);
y1 = ya1 - R*ya1*C2*sin(theta);
dya1 = sqrt((2*g*(C3*R+C4*ya1))./(A+B+(C./(D-ya1).^2).*(1+ya1./(D-ya1)).^2));
dy1  = dya1 - R*dya1*C2*sin(theta);

%% 
% Graphische Ausgabe

lgdstr = ["Näherung 1 ";"Näherung 2 ";"exakte Lsg.";"Schiene"];
           
% Vergleich Phasenraum z(x)
figure();
subplot(2,1,1)
hold on
p1(1) = plot(y,dy,'Color',Colors(3,:), 'LineWidth',2);
p1(2) = plot(y,dy,'Color',Colors(4,:), 'LineWidth',2);
p1(3) = plot(yLGLex,dyLGLex,'Color',Colors(2,:), 'LineWidth',2);
axis([0 1.2*max(yLGLex) 0 1.2*max(dyLGLex)]);
grid on
ylabel('dy in cm/s','FontSize',14)
xlabel('y in cm','FontSize',14)
h=title('Phasenraum Vergleich');
legend(p1(1:3),lgdstr(1:3,:),'location','eastoutside');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(2,1,2)
hold on
p1(1) = plot(z,dz,'Color',Colors(3,:), 'LineWidth',2);
p1(2) = plot(z,dz,'Color',Colors(4,:), 'LineWidth',2);
p1(3) = plot(zLGLex,dzLGLex,'Color',Colors(2,:), 'LineWidth',2);
axis([0.9*min(zLGLex) 1.1*max(zLGLex) 1.2*min(dzLGLex) 0]);
grid on
ylabel('dz in cm/s','FontSize',14)
xlabel('z in cm','FontSize',14)
h=title('Phasenraum Vergleich');
legend(p1(1:3),lgdstr(1:3,:),'location','eastoutside');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);



% Vergleich Zeitentwicklung
fig= figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);

yyaxis left;
hold on
p1(1)=plot(t,y,'Color',Colors(3,:), 'LineWidth',2,'LineStyle', Style(3));
p1(2)=plot(tLGLex,yLGLex,'Color',Colors(3,:), 'LineWidth',2,'LineStyle', Style(2));
p1(3)=plot(tLGLex,yLGLex,'Color',Colors(3,:), 'LineWidth',2,'LineStyle', Style(1));
grid on
axis([0 max(t) 0 1.1*max(y)]);
ylabel('y in cm','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung');
legend(p1(1:3),lgdstr(1:3,:),'location','southeast');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

yyaxis right;
hold on
p2(1)=plot(t,z,'Color',Colors(2,:), 'LineWidth',2, 'LineStyle', Style(3));
p2(2)=plot(tLGLex,zLGLex,'Color',Colors(2,:), 'LineWidth',2,'LineStyle', Style(2));
p2(3)=plot(tLGLex,zLGLex,'Color',Colors(2,:), 'LineWidth',2,'LineStyle', Style(1));
grid on
axis([0 1.1*t(kend) 0.9*min(z) 1.1*max(z)]);
legend(p2(1:3),lgdstr(1:3,:),'location','southeast');
legend box off
ylabel('z in cm','FontSize',14)
xlabel('t in s','FontSize',14)
set(gca,'FontSize',16);


figure()
hold on
p(1) = plot(y,z,'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(2));
p(2) = plot(yLGLex,zLGLex,'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(3));
p(3) = plot(yLGLex,zLGLex,'Color',Colors(2,:), 'LineWidth',2);
p(4) = plot(yr,zr,'Color',Colors(4,:), 'LineWidth',3);
PlotCircle1 (y(1),z(1),0.1*R,Colors(3,:),1, Style(3))
PlotCircle1 (y(1),z(1),R, Colors(3,:),2, Style(1))
PlotCircle1 (y(kend),z(kend),0.1*R,Colors(3,:),2,Style(1))
PlotCircle1 (y(kend),z(kend),R,Colors(3,:),1,Style(3))
axis square
axis equal
axis([-R max(y)+R -R 2*R]);
grid on
ylabel('z in cm','FontSize',14)
xlabel('x in cm','FontSize',14)
h=title('Numerische Lösung Lagrange Gl.');
legend(p,lgdstr,'location','southoutside');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

% Lagrangegleichung
function dY = dgl_DoppelKegel(t,Y,P1)
    % Y(1)- Position x(t), 
    % Y(2)- Geschwindigkeit dx(t)
    ya  = Y(1);
    dya = Y(2);
    m     = P1.m;
    R     = P1.R;
    H     = P1.H;
    g     = P1.g;
    J     = P1.J;
    A     = P1.A;
    B     = P1.B;
    C     = P1.C;
    D     = P1.D;
    C1    = P1.C1;
    C2    = P1.C2;
    C3    = P1.C3;
    C4    = P1.C4;
    C5    = P1.C5;
    phi   = P1.phi;

    dY = [dya;...
      -(2*C*D^2*dya^2 + C5*D^5*g - C5*g*ya^5 + 5*C5*D*g*ya^4 - ...
      5*C5*D^4*g*ya - 10*C5*D^2*g*ya^3 + 10*C5*D^3*g*ya^2)/...
      ((D - ya)*(A*D^4 + B*D^4 + C*D^2 + A*ya^4 + B*ya^4 + ...
      6*A*D^2*ya^2 + 6*B*D^2*ya^2 - 4*A*D*ya^3 - 4*A*D^3*ya -...
      4*B*D*ya^3 - 4*B*D^3*ya))];
end


function PlotCircle1 (xM,yM,rC,Col,LW, LS)
    u  = linspace(0,360,360);
    nx = zeros(360);
    ny = zeros(360);
    nx = rC*cosd(u)+xM;
    ny = rC*sind(u)+yM;
    plot(nx,ny,'LineWidth',LW,'Color',Col,'LineStyle',LS);
end

function [value,isterminal,direction] = MyEvent(t,Y,P1)
    value = (0.99*P1.H-Y(1)*tan(P1.phi)); 
    % erkenne max x, Abbruch, da Ungenauigkeit zu groß wird
    isterminal = 1;      % stop Integration
    direction = 0;       % negative Richtung
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


