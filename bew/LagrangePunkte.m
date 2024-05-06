% -------------------------------------------------------------------------
% LagrangePunkte.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Hill-Kurve und Lagrange-Punkte
% 
% Beispielberechnungen zur Hill-Kurve und den Lagrange-Punkten
% a) hypothethisches Massenverhältnis 1:10
% b) Jupiter-Sonne
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Parameter Fall m2/m1 = 40

m1    = 1;                                   % Masse m1
m2    =   15*m1;                             % Masse m2
a     = 1;                                   % Abstand Sonne Erde

G     = 1; 
NPkt  = 1001;                                % Anzahl Punkte
q     = m1/(m1+m2);

x1    = (q-1)*a;
x2    =  q*a;

d     = sqrt((x2 - x1)^2);
x3    = linspace(-1.5*a,1.5*a,NPkt);         % Wertebereich
y3    = linspace(-1.5*a,1.5*a,NPkt);         % Wertebereich

[X3,Y3] = meshgrid(x3,y3);

d1 = sqrt((X3-(q-1)*d).^2+Y3.^2);
d2 = sqrt((X3-q*d).^2+Y3.^2);
% Effektives Potential auf m3=1 und G=1 gesetzt:
Ueff=-m1./d1-m2./d2-(m1+m2)*(X3.^2+Y3.^2)/2/d^3;  

% Bestimmung Lagrange-Punkte 1-5 aus Ueff
TempUeff    = Ueff;
[Max3,Pos3] = max(TempUeff(501,:));
for k= 500:NPkt
  TempUeff(501,k)= NaN;
end
[Max2,Pos2] = max(TempUeff(501,:),[],'omitnan');
for k= 1:floor((x1+1.5*a)*NPkt/3/a)
  TempUeff(501,k)= NaN;
end
[Max1,Pos1] = max(TempUeff(501,:),[],'omitnan');
xL(1) = x3(Pos1);
xL(2) = x3(Pos2);
xL(3) = x3(Pos3);
xL(4) = (2*q-a)/2;
xL(5) = (2*q-a)/2;
yL(1) = 0;
yL(2) = 0;
yL(3) = 0;
yL(4) = a*sqrt(3)/2;
yL(5) = -yL(4);


%%
% Graphik
levels = linspace(Max3*0.95,Max1*1.05,20);
figure();
title('Hill-Kurve');
contour(X3,Y3,-Ueff,-levels);
hold on
PlotCircle (x2,0,0.1,Colors(10,:),1);
PlotCircle (x1,0,0.05,Colors(5,:),1);
line([xL(2) xL(3)],[0 0],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x1 xL(4)],[0 yL(4)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x2 xL(4)],[0 yL(4)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x1 xL(5)],[0 yL(5)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x2 xL(5)],[0 yL(5)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([xL(4) xL(5)],[yL(4) yL(5)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
for k = 1:5
    plot(xL(k),yL(k),'+','color',Colors(4,:),'LineWidth',2,'MarkerSize',5);
end
axis square;
grid on;
ylabel('\it{y}\rm{_3}','FontSize',14);
xlabel('\it{x}\rm{_3}','FontSize',14);
h=title('Lagrange-Punkte \it{m}\rm_2 / \it{m}\rm_1 = 40');
h.FontSize = 14;
h.FontWeight = 'normal';
set(gca,'Fontsize', 16);

%%
%Parameter Fall m2/m1 = 1000 (Jupiter)

m1    = 1;                               %Masse m1
m2    = 1000*m1;                         %Masse m2
a     = 1;                               %Abstand Sonne Erde

G     = 1; 
NPkt  = 1001;                            %Anzahl Punkte
q     = m1/(m1+m2);

x1    = (q-1)*a;
x2    =  q*a;

d     = sqrt((x2 - x1)^2);
x3    = linspace(-1.5*a,1.5*a,NPkt);         %Wertebereich
y3    = linspace(-1.5*a,1.5*a,NPkt);         %Wertebereich

[X3,Y3] = meshgrid(x3,y3);

d1 = sqrt((X3-(q-1)*d).^2+Y3.^2);
d2 = sqrt((X3-q*d).^2+Y3.^2);
% Effektives Potential auf m3=1 und G=1 gesetzt:
Ueff=-m1./d1-m2./d2-(m1+m2)*(X3.^2+Y3.^2)/2/d^3;  

% Bestimmung Lagrange-Punkte 1-5 aus Ueff
TempUeff    = Ueff;
[Max3,Pos3] = max(TempUeff(501,:));
for k= 500:NPkt
  TempUeff(501,k)= NaN;
end
[Max2,Pos2] = max(TempUeff(501,:),[],'omitnan');
for k= 1:floor((x1+1.5*a)*NPkt/3/a)
  TempUeff(501,k)= NaN;
end
[Max1,Pos1] = max(TempUeff(501,:),[],'omitnan');
xL(1) = x3(Pos1);
xL(2) = x3(Pos2);
xL(3) = x3(Pos3);
xL(4) = (2*q-a)/2;
xL(5) = (2*q-a)/2;
yL(1) = 0;
yL(2) = 0;
yL(3) = 0;
yL(4) = a*sqrt(3)/2;
yL(5) = -yL(4);


%%
% Graphik
levels = linspace(Max3*0.999,Max1*1.001,20);
figure();
contour(X3,Y3,-Ueff,-levels);
hold on
PlotCircle (x2,0,0.05,Colors(10,:),1);
PlotCircle (x1,0,0.02,Colors(5,:),1);
line([xL(2) xL(3)],[0 0],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x1 xL(4)],[0 yL(4)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x2 xL(4)],[0 yL(4)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x1 xL(5)],[0 yL(5)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([x2 xL(5)],[0 yL(5)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
line([xL(4) xL(5)],[yL(4) yL(5)],'color',Colors(4,:),...
'LineWidth',2,'LineStyle',Style(3));
for k = 1:5
    plot(xL(k),yL(k),'+','color',Colors(4,:),'LineWidth',2,'MarkerSize',5);
end
axis square;
grid on;
ylabel('\it{y}\rm{_3}','FontSize',14);
xlabel('\it{x}\rm{_3}','FontSize',14);
h=title('Lagrange-Punkte Sonne-Jupiter');
h.FontSize = 14;
h.FontWeight = 'normal';
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

