% -------------------------------------------------------------------------
% Wandsonnenuhr.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Lineatur einer Wandsonnenuhr
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":","-", "-.", ":", "--", ":"];

%% Numerische Lösung
% Parameter

phi     =  deg2rad(48);     % Geographische Breite
R       =  1;
sigmaA   = deg2rad(linspace(-90,90,7));
sigmaV   = atan(cos(phi)*tan(sigmaA));
hours    = (rad2deg(sigmaA)+90)/15+6;

figure();
subplot(1,2,1)
hold on;
delta   =  deg2rad(-22.5);     % Wandverdrehung 
sigmaVs  = atan(cos(phi)./(cos(delta).*(1./tan(sigmaA))+sin(phi)*sin(delta)));
str1 = strcat('Ausrichtung ONO-WSW : (', string(num2str(rad2deg(delta),3)));
str1 = strcat(str1,'°)');
ttl=title(str1);
for k = 1:length(sigmaA) 
    x = sin(sigmaV(k));
    y = -cos(sigmaV(k));
    q(k)=line([0 x],[0 y],...
         'color',Colors(4,:),'LineWidth',1,'LineStyle',Style(3));
    text(x-0.02,y-0.02,num2str(hours(k),2),'color',Colors(4,:));
end
for k = 1:length(sigmaA) 
    x = sin(sigmaVs(k));
    y = -cos(sigmaVs(k));
    if sigmaVs(k)*sigmaV(k) >= 0 
        p(k)=line([0 x],[0 y],...
         'color',Colors(2,:),'LineWidth',2,'LineStyle',Style(1));
        text(x-0.02,y-0.02,num2str(hours(k),2),'color',Colors(2,:));
    end
end
grid on
axis equal
axis equal
ylim([-1.2 0]);
xlim([-1.2 1.2]);
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);


subplot(1,2,2)
delta   =  deg2rad(+22.5);     % Wandverdrehung 
sigmaVs  = atan(cos(phi)./(cos(delta).*(1./tan(sigmaA))+sin(phi)*sin(delta)));
hold on;
str1 = strcat('Ausrichtung OSO-WNW : (', string(num2str(rad2deg(delta),3)));
str1 = strcat(str1,'°)');
ttl= title(str1);

for k = 1:length(sigmaA) 
    x = sin(sigmaV(k));
    y = -cos(sigmaV(k));
    q(k)=line([0 x],[0 y],...
         'color',Colors(4,:),'LineWidth',1,'LineStyle',Style(3));
    text(x-0.02,y-0.02,num2str(hours(k),2),'color',Colors(4,:));
end
for k = 1:length(sigmaA) 
    x = sin(sigmaVs(k));
    y = -cos(sigmaVs(k));
    if sigmaVs(k)*sigmaV(k) >= 0 
        p(k)=line([0 x],[0 y],...
         'color',Colors(3,:),'LineWidth',2,'LineStyle',Style(1));
        text(x-0.02,y-0.02,num2str(hours(k),2),'color',Colors(3,:));
    end
end
grid on
axis equal
axis equal
ylim([-1.2 0]);
xlim([-1.2 1.2]);
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------