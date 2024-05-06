% -------------------------------------------------------------------------
% ElastischerStoss.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet den Streuwinkel fuer die Masse m_1 im Laborsystem 
% fuer verschiedene % Massenverhaeltnisse als Funktion des Streuwinkels
% im Schwerpunktsystem. (Die Masse m_2 ruht vor dem Stoss.)
% 
% Programm berechnet den Energieuebertrag als Funktion des Streuwinkels
% im Laborsystem.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;


% Anfangsbedingungen, Parameter
mratio = [0.01, 0.1, 0.5, 0.75, 1, 1.1, 2, 5 8];
Stil   = ["--", "--", "--", "--", '-', ':', ':', ':', ':'];
thetaS = linspace(0,deg2rad(180),1000);

for k=1:length(mratio)
    theta1(k,:) = atan2d(sin(thetaS),(cos(thetaS)+mratio(k)));
end

figure();
for k=1:length(mratio)
    plot(rad2deg(thetaS),theta1(k,:),'Color',Colors(k,:),...
        'LineStyle',Stil(k),'LineWidth',1);
    hold on;
end    
xlabel('Streuwinkel {\theta}_s in °','FontSize',12);
ylabel('Streuwinkel {\theta}_1 in °','FontSize',12);
legend(strcat('m_1/m_2= ', num2str(mratio(:),'%4.2f')),'location','northwest');
legend box off
axis square;
grid on
Ticks = 0:30:180;
set(gca, 'XTickMode', 'manual', 'XTick', Ticks, 'xlim', [0,180]);
set(gca, 'YTickMode', 'manual', 'YTick', Ticks, 'ylim', [0,180]);
set(gca,'Fontsize', 16);

for k=1:length(mratio)
    fac = mratio(k)/(1+mratio(k)); 
    ET(k,:) = 2*fac*(1-fac)*cos(thetaS);
end

figure();
for k=1:length(mratio)
  subplot(1,2,1);
  plot(rad2deg(thetaS),ET(k,:)*100,'Color',Colors(k,:),...
     'LineStyle',Stil(k),'LineWidth',1);
  hold on;
end    
xlabel('Streuwinkel {\theta}_s in °','FontSize',12);
ylabel('Energieuebertrag in %','FontSize',12);
legend(strcat('m_1/m_2= ',num2str(mratio(:),'%4.2f')),'location','northeast');
legend box off
grid on
set(gca, 'XTickMode', 'manual', 'XTick', Ticks, 'xlim', [0,180]);
set(gca,'Fontsize', 16);

for k=1:length(mratio)
  subplot(1,2,2);
  plot(theta1(k,:),ET(k,:),'Color',Colors(k,:),...
     'LineStyle',Stil(k),'LineWidth',1);
  hold on;
end    
xlabel('Streuwinkel {\theta}_1 in °','FontSize',12);
ylabel('Energieuebertrag in %','FontSize',12);
legend(strcat('m_1/m_2= ',num2str(mratio(:),'%4.2f')),'location','northeast');
legend box off
grid on
set(gca, 'XTickMode', 'manual', 'XTick', Ticks, 'xlim', [0,180]);
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




