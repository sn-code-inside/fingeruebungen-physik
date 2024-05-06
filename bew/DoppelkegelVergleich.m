% -------------------------------------------------------------------------
% DoppelkegelVergleich.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärts rollender Doppelkegel
% 
% Programm zeigt korrekte und inkorrekte bzw. nur näherungsweise 
% gültige Lösungen für den aufwärts rollenden Doppelkegel  
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


%%
% Import numerische Berechnung aus Doppelkegel02.m (LGL)
QLGL = ImportFile01('Data/DKkorrekt.txt');
tLGLex  =  QLGL.t;
yLGLex  =  QLGL.y;
zLGLex  =  QLGL.z;
dyLGLex =  QLGL.dy;
dzLGLex =  QLGL.dz;

% Import numerische Berechnung aus Doppelkegel05.m (LGL) Näherung1
QLGL1 = ImportFile01('Data/DKApprox01.txt');
t1  =  QLGL1.t;
y1  =  QLGL1.y;
z1  =  QLGL1.z;
dy1 =  QLGL1.dy;
dz1 =  QLGL1.dz;

% Import numerische Berechnung aus Doppelkegel06.m (LGL) Näherung2
QLGL2 = ImportFile01('Data/DKApprox02.txt');
t2  =  QLGL2.t;
y2  =  QLGL2.y;
z2  =  QLGL2.z;
dy2 =  QLGL2.dy;
dz2 =  QLGL2.dz;


% Schiene
theta   = deg2rad(06.5);   % Anstieg Schiene
R       = 3;               % Radius Kegel
yr = linspace(0,25);
zr = yr*tan(theta);


%% 
% Graphische Ausgabe

lgdstr = ["Näherung 1 ";"Näherung 2 ";"exakte Lsg.";"Schiene"];
           
% Vergleich Phasenraum z(x)
figure();
subplot(2,1,1)
hold on
p1(1) = plot(y1,dy1,'Color',Colors(3,:), 'LineWidth',2);
p1(2) = plot(y2,dy2,'Color',Colors(4,:), 'LineWidth',2);
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
p1(1) = plot(z1,dz1,'Color',Colors(3,:), 'LineWidth',2);
p1(2) = plot(z2,dz2,'Color',Colors(4,:), 'LineWidth',2);
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
p1(1)=plot(t1,y1,'Color',Colors(3,:), 'LineWidth',2,'LineStyle', Style(1));
p1(2)=plot(t2,y2,'Color',Colors(4,:), 'LineWidth',2,'LineStyle', Style(1));
p1(3)=plot(tLGLex,yLGLex,'Color',Colors(2,:), 'LineWidth',2,'LineStyle', Style(1));
grid on
axis([0 max(tLGLex) 0 1.1*max(yLGLex)]);
ylabel('y in cm','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Zeitentwicklung');
legend(p1(1:3),lgdstr(1:3,:),'location','southeast');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

yyaxis right;
hold on
p2(1)=plot(t1,z1,'Color',Colors(3,:), 'LineWidth',2, 'LineStyle', Style(2));
p2(2)=plot(t2,z2,'Color',Colors(4,:), 'LineWidth',2,'LineStyle', Style(2));
p2(3)=plot(tLGLex,zLGLex,'Color',Colors(2,:), 'LineWidth',2,'LineStyle', Style(2));
grid on
axis([0 1.1*max(tLGLex) 0.9*min(zLGLex) 1.1*max(zLGLex)]);
legend(p2(1:3),lgdstr(1:3,:),'location','southeast');
legend box off
ylabel('z in cm','FontSize',14)
xlabel('t in s','FontSize',14)
set(gca,'FontSize',16);

% Geometrie
figure()
hold on
p(1) = plot(y1,z1,'Color',Colors(3,:), 'LineWidth',2,'LineStyle',Style(2));
p(2) = plot(y2,z2,'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(3));
p(3) = plot(yLGLex,zLGLex,'Color',Colors(2,:), 'LineWidth',2);
p(4) = plot(yr,zr,'Color',Colors(4,:), 'LineWidth',3);
PlotCircle1 (yLGLex(1),zLGLex(1),0.1*R,Colors(3,:),1, Style(3))
PlotCircle1 (yLGLex(1),zLGLex(1),R, Colors(3,:), 2, Style(1))
PlotCircle1 (yLGLex(end),zLGLex(end),0.1*R,Colors(3,:),2,Style(1))
PlotCircle1 (yLGLex(end),zLGLex(end),R,Colors(3,:),1,Style(3))
axis square
axis equal
axis([-R max(yr)+R -R 2*R]);
grid on
ylabel('z in cm','FontSize',14)
xlabel('y in cm','FontSize',14)
h=title('Numerische Lösung Lagrange Gl.');
legend(p,lgdstr,'location','southoutside');
legend box off
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% 
% Funktionen

function PlotCircle1 (xM,yM,rC,Col,LW, LS)
    u = linspace(0,360,360);
    nx = zeros(360);
    ny = zeros(360);
    nx = rC*cosd(u)+xM;
    ny = rC*sind(u)+yM;
    plot(nx,ny,'LineWidth',LW,'Color',Col,'LineStyle',LS);
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


