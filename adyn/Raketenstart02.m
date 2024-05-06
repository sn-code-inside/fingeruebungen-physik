% -------------------------------------------------------------------------
% Raketenstart02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Nutzlast-Antriebsvermögen-Relation als Funktion des
% Strukturmassen-Verhältnisses
%
%--------------------------------------------------------------------------


%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;

%%
% Initialisierung
% alle Daten in kg, m und s 


sigma  = [0.05, 0.1, 0.15, 0.2, 0.25];
vG1    = 3000;
vG2    = 4000;
muL     = linspace(0.0001,1,10000);

for k=1:length(sigma)
   ratio_V(k,:) = log(1./(sigma(k)+muL(:))); 
   lgdstr(k,:) = num2str(sigma(k),'%3.2f');
end

%%
% Graphische Ausgabe

figure()
% Plot Nutzlastverhältnis über Antriebsvermögen
for k=1:length(sigma)
    semilogy(ratio_V(k,:)*vG2, muL(:),'Color', Colors(k,:),'LineWidth',2);
    hold on
end
grid on;
xlim([0, 3*vG2]);
ylabel('Nutzlastverhältnis','FontSize',14)
xlabel('Antriebsvermögen Delta v','FontSize',14)
h2 = legend(lgdstr,'location','northeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
set(gca,'FontSize',16);

