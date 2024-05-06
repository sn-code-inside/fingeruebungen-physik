% -------------------------------------------------------------------------
% Doppelkegel02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärtsrollender Doppelkegel
% 
% Programm berechnet die Lösung der
% Lagrange-Gleichungen des aufwärtsrollenden Doppelkegels  
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerische Lösung

% Parameter

H       = 6.5;              % Länge Schaukel cm
R       = 3.0;              % Länge Oberkörper cm 
phi     = deg2rad(15.3);    % Öffnungswinkel Schienen
theta   = deg2rad(06.5);    % Anstiegswinkel Schienen
psi     = atan(R/H);        % Öffnungswinkel Kegel
M       = 122.5;            % Masse Doppelkegel in g
J       = 3*M*R^2/10;       % Trägheitsmoment Doppelkegel in g cm²
g       = 981;              % g in cm/s

alpha   = asin(tan(psi)*tan(phi));  % Richtungswinkel alpha
q0      = 0;                % Anfangsgposition gen. KO des Doppelkegels
dq0     = 0;                % Anfangsgeschwindigkeit gen. KO des Doppelkegels
ys0     = H*tan(psi)*sin(alpha);
                            % Anfangsgposition ys im KOS Sigma
                         
fprintf('\n ');
fprintf('\n phi   = %4.2f°', rad2deg(phi));
fprintf('\n psi   = %4.2f°', rad2deg(psi));
fprintf('\n theta = % 4.2f°', rad2deg(theta));
fprintf('\n alpha = % 4.2f°', rad2deg(alpha));
fprintf('\n ');
fprintf('\n H  = %8.2f cm', H);
fprintf('\n R  = %8.2f cm', R);
fprintf('\n M  = %8.2f g', M);
fprintf('\n J  = %8.2f gcm²', J);
fprintf('\n ');
fprintf('\n ys0= %8.2f cm', ys0);
fprintf('\n g  = %8.2f cm/s', g);
fprintf('\n ');
fprintf('\n Rollbedingung erfüllt ??? ');
b1 = tan(theta);
b2 = tan(phi)*tan(psi)/sqrt(1-(tan(phi)*tan(psi))^2);
if b1 < b2
    fprintf('\n %s  < %s ! Rollbedingung erüllt!!!',num2str(b1,4), ...
        num2str(b2,4));
    fprintf('\n ');
else
    fprintf('\n %s  > %s ! Rollbedingung nicht(!) erüllt! ',...
        num2str(b1,4), num2str(b2,4));
    fprintf('\n ');
end
fprintf('\n ');

%% Berechnungen

% Anfangswerte
AB=[q0;dq0]; % AB für ode45
tmax  = 25;  % maximale Berechnungszeit in s

% Parameterset für ODE45
P1.H =H;
P1.R =H;
P1.alpha = alpha;
P1.psi   = psi;
P1.phi   = phi;
P1.C1 = 3*R*R/10/(tan(alpha))^2;
P1.C2 = -g*sin(alpha-theta);
P1.C3 = H*tan(psi)/tan(alpha);

% Numerische Lösung Lagrangegleichung
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6,'events',@MyEvent);
% Numerische Lösung LGL
[t, Y, TE, YE, IE] = ode45(@dgl_DoppelKegel,[0.0,tmax], AB,opt,P1); 
q  = Y(:,1);
dq = Y(:,2);
ys = ys0 + Y(:,1)*cos(alpha);
dys = Y(:,2)*cos(alpha);
zs = -ys*tan(alpha) + H*tan(psi)*(cos(alpha)-tan(alpha)*sin(alpha));
dzs = -dys*tan(alpha);

kend = length(t);
tend = TE;

% Transformation ins KOS Sigma
y  = ys*cos(theta) - zs*sin(theta);
z  = ys*sin(theta) + zs*cos(theta);
dy  = dys*cos(theta) - dzs*sin(theta);
dz  = dys*sin(theta) + dzs*cos(theta);

% Schienen
yr = linspace(0,25);
zr = yr*tan(theta);
% Vektor a(q)
a  = H*tan(psi)-q*tan(alpha);

fileID = fopen('Data/DKkorrekt.txt','w');
fprintf(fileID,'%s %s %s %s %s \r\n','t','y','dy','z', 'dz');
for k=1:length(t) 
    fprintf(fileID,'%10.4f %10.4f %10.4f %10.4f %10.4f\r\n',...
                    t(k),y(k),dy(k),z(k),dz(k));
end
fclose(fileID);

%% Graphische Ausgabe

% Lösung Lagrangegleichung z(x)
fig1=figure();
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
h=title('Numerische Lösung Lagrange Gl.');
legend('Schwerpunkt S','Schiene','location','southeast');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);


% Phasenraum 
fig2= figure();
set(fig2,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
yyaxis left;
hold on
plot(t,y,'Color',Colors(3,:), 'LineWidth',2);
grid on
axis([0 1.1*tend 0 1.1*max(y)]);
ylabel('y in cm','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

yyaxis right;
hold on
plot(t,z,'Color',Colors(2,:), 'LineWidth',2, 'LineStyle', Style(1));
grid on
axis([0 1.1*tend 0.9*min(z) 1.1*max(z)]);
ylabel('z in cm','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung y, z');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

fig3=figure();
set(fig3,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
yyaxis left;
hold on
plot(t,q,'Color',Colors(3,:), 'LineWidth',2);
grid on
axis([0 1.1*t(kend) 0 1.1*max(q)]);
ylabel('q in cm','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung q, dq');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

yyaxis right;
hold on
plot(t,dq,'Color',Colors(2,:), 'LineWidth',2, 'LineStyle', Style(1));
grid on
axis([0 1.1*t(kend) 0 1.2*max(dq)]);
ylabel('dq in cm/s','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung q, dq/dt');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

% Lagrangegleichung
function dY = dgl_DoppelKegel(t,Y,P1)
    % Y(1)- Position q(t), 
    % Y(2)- Geschwindigkeit dq(t)
    Q     = 1-P1.C2*Y(1);
    dY    = [Y(2);...    
         -(P1.C2*P1.C3^3 + P1.C1*Y(2)^2 - P1.C2*Y(1)^3 + ...
           3*P1.C2*P1.C3*Y(1)^2 - 3*P1.C2*P1.C3^2*Y(1))/...
         ((P1.C3 - Y(1))*(P1.C1 - 2*P1.C3*Y(1) + P1.C3^2 + Y(1)^2))];
end

function PlotCircle1 (xM,yM,rC,Col,LW, LS)
    u=linspace(0,360,360);
    nx=zeros(360);
    ny=zeros(360);
    nx= rC*cosd(u)+xM;
    ny= rC*sind(u)+yM;
    plot(nx,ny,'LineWidth',LW,'Color',Col,'LineStyle',LS);
end

function [value,isterminal,direction] = MyEvent(t,Y,P1)
    value = (0.99*P1.H*tan(P1.psi)/tan(P1.alpha)-Y(1)); 
    % erkenne max q, Abbruch, da Ungenauigkeit zu groß wird
    isterminal = 1;      % stop Integration
    direction = 0;       % negative Richtung
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

