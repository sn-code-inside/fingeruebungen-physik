% -------------------------------------------------------------------------
% Doppelkegel05.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärts rollender Doppelkegel
% 
% Programm berechnet inkorrekte nur näherungsweise 
% gültige Lagrange-Funktion (Näherung 1) und 
% Lagrange-Gleichungen des aufwärts rollenden Doppelkegels  
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

syms m J R  y dy A B C1 C2 C3 C4 g 'real'


%% Symbolische Lösung Lagrange-Funktion

% Verallgemeinerte Koordinaten und Ableitungen
Q  = [y];
dQ = [dy];


% kinetische Energie
T_kin = 1/2*m*dy.*dy.*...
    (1 + C2^2 + (C4./(R-C1.*y)^2).*(1 + C1*y./(R-C1.*y)).^2);

% potentielle Energien
U_pot =  m*g*(R+y*C2);


% Lagrange-Funktion, Energie
L = T_kin - U_pot;
E = T_kin + U_pot;

%% Ausgabe der Ergebnisse
fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf(1, 'E_0 = T + U = %s \n', E);
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
C2    = tan(theta)-C1;
C3   = cos(theta);
C4    = J/m/C3^2;

y0 = 0;
dy0= 0;


%% Berechnungen ODE45

% Anfangswerte
AB=[y0;dy0]; % AB für ode45

P1.m     = m;
P1.g     = g;
P1.J     = J;
P1.H     = H;
P1.R     = R;
P1.C1    = C1;
P1.C2    = C2;
P1.C3    = C3;
P1.C4    = C4;
P1.A     = A;
P1.B     = B;
P1.phi   = phi;


% Numerische Lösung Lagrangegleichung
opt=odeset('AbsTol',1.e-10,'RelTol',1.e-10,'events',@MyEvent);
% Numerische Lösung LGL
[t, Y, TE, YE, IE] = ode45(@dgl_DoppelKegel,[0.0,tmax],AB,opt,P1); 
y  = Y(:,1);
dy = Y(:,2); 
z = R + y*C2;
dz = dy*C2;
kend = length(y);
yr = linspace(0,25);
zr = yr*tan(theta);
E0 = 0;
y1 = linspace(0,max(y),100);
dy1= sqrt((-2*g*C2*y1)./(1+C2^2+C4./(R-C1.*y1).^2).*...
          (1+C1*y1./(R-C1*y1)).^2);
      
fileID = fopen('Data/DKApprox01.txt','w');
fprintf(fileID,'%s %s %s %s %s \r\n','t','y','dy','z', 'dz');
for k=1:length(t) 
    fprintf(fileID,'%10.4f %10.4f %10.4f %10.4f %10.4f\r\n',...
                    t(k),y(k),dy(k),z(k),dz(k));
end
fclose(fileID);


%% 
% Graphische Ausgabe

% Phasenraum 
figure();
subplot(2,1,1)
hold on
p1(1) = plot(y,dy,'Color',Colors(3,:), 'LineWidth',2);
axis([0 1.2*max(y) 0 1.2*max(dy)]);
grid on
ylabel('dy in cm/s','FontSize',14)
xlabel('y in cm','FontSize',14)
h=title('Phasenraum y');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(2,1,2)
hold on
p1(1) = plot(z,dz,'Color',Colors(2,:), 'LineWidth',2);
axis([0.9*min(z) 1.1*max(z) 1.2*min(dz) 0]);
grid on
ylabel('dz in cm/s','FontSize',14)
xlabel('z in cm','FontSize',14)
h=title('Phasenraum z');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);




% Zeitentwicklung 
fig= figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);

yyaxis left;
hold on
plot(t,y,'Color',Colors(3,:), 'LineWidth',2);
grid on
axis([0 max(t) 0 1.1*max(y)]);
ylabel('y in cm','FontSize',14)
xlabel('t in s','FontSize',14)
set(gca,'FontSize',16);

yyaxis right;
hold on
plot(t,z,'Color',Colors(2,:), 'LineWidth',2, 'LineStyle', Style(1));
grid on
axis([0 1.1*t(kend) 0.9*min(z) 1.1*max(z)]);
ylabel('z in cm','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);



% Geometrie z(x)
figure();
hold on
plot(y,z,'Color',Colors(3,:), 'LineWidth',2);
plot(yr,zr,'Color',Colors(4,:), 'LineWidth',3);
PlotCircle1 (y(1),z(1),0.1*R,Colors(3,:),1, Style(3))
PlotCircle1 (y(1),z(1),R, Colors(3,:),2, Style(1))
PlotCircle1 (y(kend),z(kend),0.1*R,Colors(3,:),2,Style(1))
PlotCircle1 (y(kend),z(kend),R,Colors(3,:),1,Style(3))
axis square
axis equal
axis([-R max(y)+R -R 2*R]);
grid on
ylabel('z in cm','FontSize',14)
xlabel('y in cm','FontSize',14)
h1=title('Numerische Lösung Lagrange Gl.');
set(h1,'FontSize',12,'FontWeight','normal'); 
legend('Näherung 1','Schiene','location','southeast');
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




%% Funktionen

% Lagrangegleichung
function dY = dgl_DoppelKegel(t,Y,P1)
    % Y(1)- Position x(t), 
    % Y(2)- Geschwindigkeit dx(t)
    y  = Y(1);
    dy = Y(2);
    m     = P1.m;
    R     = P1.R;
    H     = P1.H;
    g     = P1.g;
    J     = P1.J;
    A     = P1.A;
    B     = P1.B;
    C1    = P1.C1;
    C2    = P1.C2;
    C3    = P1.C3;
    C4    = P1.C4;
    phi   = P1.phi;

    dY = [dy;...
      -(C2*R^5*g + 2*C1*C4*R^2*dy^2 - C1^5*C2*g*y^5 + ...
      10*C1^2*C2*R^3*g*y^2 - 10*C1^3*C2*R^2*g*y^3 - 5*C1*C2*R^4*g*y +...
      5*C1^4*C2*R*g*y^4)/((R - C1*y)*(C4*R^2 + R^4 + C2^2*R^4 + ...
      C1^4*y^4 - 4*C1^3*R*y^3 + C1^4*C2^2*y^4 + 6*C1^2*R^2*y^2 -...
      4*C1*R^3*y + 6*C1^2*C2^2*R^2*y^2 -...
      4*C1*C2^2*R^3*y - 4*C1^3*C2^2*R*y^3))];
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
    value = (25.0-Y(1)); 
    % erkenne max x, Abbruch, da Ungenauigkeit zu groß wird
    isterminal = 1;      % stop Integration
    direction = 0;       % negative Richtung
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


