% -------------------------------------------------------------------------
% UmlaufErdeMond.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Bewegung von Punkten der Erde bei Umlauf um Mond.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
xS = 0;
yS = 0;

RS = 4675;
RE = 6370;
NoDays = 36;

days = linspace(0,NoDays,NoDays+1);
%Mittelpunkt
yM = RS*cosd(days*360/NoDays)/RE;
xM = RS*sind(-days*360/NoDays)/RE;
%Top
yT = RS*cosd(days*360/NoDays)/RE+1;
xT = RS*sind(-days*360/NoDays)/RE;
%Bottom
yB = RS*cosd(days*360/NoDays)/RE-1;
xB = RS*sind(-days*360/NoDays)/RE;

%Richtung Mond
xMoon = 2*sind((days)*360/NoDays);
yMoon = -2*cosd((days)*360/NoDays);

Delay = 0.05;
figure()
%Schwerpunkt Erde-Mond
PlotCircle (xS,yS,0.02,Colors(4,:),5);
hold on
axis equal;
grid on;
axis equal
ylim([-2.5 2.5]);
xlim([-2.5 2.5]);
xlabel('a.u.');
ylabel('a.u.');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);
pause(Delay);
%Bewegung Erdmittelpunkt
plot(xM(1),yM(1),'+','Color',Colors(3,:),'LineWidth',1);
h = PlotCircle(xM(1),yM(1),1,Colors(3,:),1);
h1 = line([xS, xMoon(1)],[yS,yMoon(1)],...
        'Color',Colors(10,:),'LineWidth',1,'LineStyle',':');
PlotCircle (xMoon(1),yMoon(1),0.02,Colors(10,:),5);
for k=2:NoDays+1
  pause(Delay);
  set(h,'Visible','off');
%   set(h1,'Visible','off');
  plot(xM(k),yM(k),'+','Color',Colors(3,:),'LineWidth',1);
  h = PlotCircle(xM(k),yM(k),1,Colors(3,:),1);
  h1 = line([xS, xMoon(k)],[yS,yMoon(k)],...
        'Color',Colors(10,:),'LineWidth',1,'LineStyle',':');
  PlotCircle (xMoon(k),yMoon(k),0.02,Colors(10,:),5);
end
%Bewegung "Top"
plot(xT(1),yT(1),'d','Color',Colors(2,:),'LineWidth',1);
h = PlotCircle(xM(1),yM(1),1,Colors(2,:),1);
for k=2:NoDays+1
  pause(Delay);
  set(h,'Visible','off');
  plot(xT(k),yT(k),'d','Color',Colors(2,:),'LineWidth',1);
  h = PlotCircle(xM(k),yM(k),1,Colors(2,:),1);
end
%Bewegung "Bottom"
plot(xB(1),yB(1),'s','Color',Colors(4,:),'LineWidth',1);
h = PlotCircle(xM(1),yM(1),1,Colors(4,:),1);
for k=2:NoDays+1
  pause(Delay);
  set(h,'Visible','off');
  plot(xB(k),yB(k),'s','Color',Colors(4,:),'LineWidth',1);
  h = PlotCircle(xM(k),yM(k),1,Colors(4,:),1);
end
 h = PlotCircle(xM(k),yM(k),1,Colors(3,:),1);
%Schwerpunkt Erde-Mond
PlotCircle (xS,yS,0.02,Colors(4,:),5);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

function h = PlotCircle (xM,yM,rC,Col,LW)
    u=linspace(0,360,360);
    nx=zeros(360);
    ny=zeros(360);
    Nmax = 1;
    nx= rC*cosd(u)/Nmax+xM;
    ny= rC*sind(u)/Nmax+yM;
    h = plot(nx,ny,'LineWidth',LW,'Color',Col);
 end




