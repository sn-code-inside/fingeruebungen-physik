% -------------------------------------------------------------------------
% AnHarmonOszillator01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Schwingungsperiode Nichtlineares Pendel
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
NPts      = 300;                            % Schritte
thetamaxD = 90;                             % Amplitude
thetamax = deg2rad(thetamaxD); 
thetastep = thetamax/NPts;
thetaA  = 0:thetastep:thetamax;             % in rad
thetaAD =rad2deg(thetaA);
k_ell  = sin(thetaA/2);
m_ell  = k_ell.^2 ;
T_exakt= 2*ellipke(m_ell)/pi;
TN(1,:)  = 1./sqrt(1-(thetaA).^2/8);        % Näherung Skript
TN(2,:)  = 1*(1+(thetaA).^2/16);            % Näherung Bernoull
TN(3,:)  = 1./sqrt(cos(thetaA/2));          % Näherung Kidd
TN(4,:)  = 1./(1-(thetaA).^2/16);           % Näherung Denmann
plot(thetaAD,T_exakt,'Color',Colors(2,:), 'LineWidth',2);
hold on
for k=1:4
    plot(thetaAD,TN(k,:),'Color',Colors(2+k,:),'LineStyle', Style(k));
end
line([0,thetamaxD],[1,1],'Color',Colors(1,:),'LineStyle', Style(2)); % Harmon. Oszillator
axis([0 thetamaxD 0.95 1.2]);
h=legend('Ellipt. Integr.', 'Näherung Skript', 'Bernoulli',...
         'Kidd','Denmann','Harmon. Oszill.', 'location','northwest' );
set(h,'FontSize',14)
legend box off
xlabel('Amplitude(°)','FontSize',14)
ylabel('Schwingungsperiode','FontSize',14)
grid on
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
