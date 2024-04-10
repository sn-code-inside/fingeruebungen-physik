% -------------------------------------------------------------------------
% RocheLimit.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Roche-Limit und Hill-Sphären für Planeten 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Import Data
% alle Daten in kg, m und s
PlanetPara = ImportPlanetParameter('PlanetenParameter.csv')

abst   = PlanetPara.aP;
Masse  = PlanetPara.Masse;
rho    = PlanetPara.rhoP;
RadiusP= PlanetPara.RP;

rhoKom = PlanetPara.rhoP(11);
AE     = PlanetPara.aP(10);
MS     = PlanetPara.Masse(10);
G      = 6.672e-11;

% Begin Rechnungen
for i=1:9
    HillSph(i)    = abst(i)*(Masse(i)/3/MS)^(1/3);
    HillSphAE(i)  = HillSph(i)/AE;
    RocheLim(i)   = 1.26*RadiusP(i)*(rho(i)/rhoKom)^(1/3);
end

LgdStr = PlanetPara.Name;


%--------------------------------------------------------------------------
% Graphische Ausgabe
figure()
hold on;
for iPlot = 1:9
  h(iPlot) = plot(iPlot,RocheLim(iPlot)/1000,'d','Color', Colors(iPlot,:),...
             'LineWidth',2);
  hold on;
end
grid on;
header2 = sprintf('Roche-Limit');
ttl = title(header2);
ttl.FontSize = 16;
ttl.FontWeight = 'normal';
xlabel('Planet Nr')
ylabel('Roche-Limit in km');

lgd=legend(h,LgdStr(1:9),'Location','bestoutside','NumColumns',1);
lgd.FontSize=16;
legend boxoff;
set(gca, 'Fontsize', 14, 'linewidth', 1);

figure()
hold on;
for iPlot = 1:9
  x(iPlot) = iPlot;
  h(iPlot) = plot(iPlot,log(HillSphAE(iPlot)),'d','Color', Colors(iPlot,:),...
             'LineWidth',2);
  hold on;
end
grid on;
header2 = sprintf('Hill-Sphäre in AE');
ttl = title(header2);
ttl.FontSize = 16;
ttl.FontWeight = 'normal';
xlabel('Planet Nr')
ylabel('Logarithmus der Hill-Sphäre in AE');

lgd=legend(h,LgdStr(1:9),'Location','bestoutside','NumColumns',1);
lgd.FontSize=16;
legend boxoff;
set(gca, 'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
