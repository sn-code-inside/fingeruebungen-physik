% -------------------------------------------------------------------------
% FoucaultPendel.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Foucault-Pendel
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%

phi  = 70;             % Breite in °
g    = 9.81;              
L    = 100;            % Pendelperiode 
Tp   = 2*pi*sqrt(L/g); % Pendellänge
omega0 = 2*pi/Tp;      % Eigenfreuenz Pendel

omegaE0  = 2*pi/86164;               % Erd-Rotationskreisfrequenz
omegaE  = omegaE0*300;               % Erd-Frequenz-Vielfaches
omegaF  = omegaE*sind(phi);          % Foucault Frequenz
omegaD  = sqrt(omegaF^2+omega0^2);   % leicht verschobene Pendelfrequenz
NPoints = 5000;       % Stützstellen
x0=1.0; y0=0.0;       % AB
tF=2*pi/omegaF;       % Foucault Präzessionsperiode
tmax=tF  ;            % Halber Umlauf 
ts=tmax/NPoints;      % Zeit-Stützstellen

t=[0:ts:tmax]; NPts=length(t); nh=round(NPts/2);

% Pendelkoordinaten
x=x0*cos(omegaF*t).*cos(omegaD*t)+y0*sin(omegaF*t).*cos(omegaD*t);
y=-x0*sin(omegaF*t).*cos(omegaD*t)+y0*cos(omegaF*t).*cos(omegaD*t);

str1=cat(2,'Breite = ',num2str(phi,3),' °');
str2=cat(2,'T_F = ',num2str(tF/60,3),' min');
str3=cat(2,'\omega_E/\omega_{E0} = ',num2str(omegaE/omegaE0,3),'');
figure()
% erster Umlauf
plot(x(1:nh),y(1:nh),'Color',Colors(2,:),'LineStyle',Style(1),...
    'LineWidth',2'); 
hold on; 
% zweiter Umlauf
plot(x(nh:NPts),y(nh:NPts),'Color',Colors(4,:),'LineStyle',Style(3),...
     'LineWidth',2'); 
title('Foucault-Pendel','FontSize',14)
ylabel('y\prime','FontSize',14); xlabel('x\prime','FontSize',14);

minx= min(x)*1.2;
miny= min(y)*1.2;
maxx= max(x)*1.2;
maxy= max(y)*1.2;
axis([minx, maxx, miny,maxy])

axis square;
grid on;
text(minx*(1-0.05),miny*(1-0.10),str1,'FontSize',10,'Color','blue');
text(minx*(1-0.05),miny*(1-0.25),str2,'FontSize',10,'Color','blue');
text(minx*(1-0.05),miny*(1-0.40),str3,'FontSize',10,'Color','blue');

% Simulation  
figure();
ylabel('y\prime','FontSize',14); xlabel('x\prime','FontSize',14);
minx = min(x);
miny = min(y);
maxx = max(x);
maxy = max(y);
axis([minx, maxx, miny,maxy])
axis square;
 for k=1:3:NPoints/2
     axis([minx, maxx, miny,maxy])
     hold on
     strtime = strjoin({' t = ', num2str(t(k), '% 6.0f s')});
     h1=text(min(x)*(1-0.05),max(y)*(1-0.15),strtime,'FontSize',10,'Color','blue');
     h1.Visible='on';
     plot(x(k),y(k),'o',...
          'MarkerSize',2,...
          'MarkerEdgeColor',Colors(2,:),...
          'MarkerFaceColor',Colors(2,:)'); %m1 @ x1,y1
     pause(.005)
     h1.Visible='off';
 end
 h1.Visible='on';
ylabel('y\prime','FontSize',14); xlabel('x\prime','FontSize',14);
axis([minx, maxx, miny,maxy])
axis square;
grid on;
text(min(x)*(1-0.05),min(y)*(1-0.10),str1,'FontSize',10,'Color','blue');
text(min(x)*(1-0.05),min(y)*(1-0.25),str2,'FontSize',10,'Color','blue');
text(min(x)*(1-0.05),min(y)*(1-0.40),str3,'FontSize',10,'Color','blue');
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------