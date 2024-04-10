% -------------------------------------------------------------------------
% PerihelKappa3.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Perihelverschiebung mit Stoerterm kappa3 ungleich Null 
% (vgl. Kapitel 2.5.4).
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
ups0 = 0;
ups=linspace(0, 4*pi, 4000);    % zwei Perioden 

%Parameter
a   = 1;                        % Radius
mu_G= 1;                        % Gravitationskonstante * Masse
ex  = 0.1;                      % Exzentrizitaet

% Berechnung des Drehimpuls
L   = sqrt(a*mu_G*(1-ex*ex));
% Berechnung der Parameter aus Formel 3.131 
p   = L*L/mu_G;
C   = 1+0.5*ex^2-ex^2.*cos(2*ups)/6;

% kappa3 = 0
r0  = p./(1+ex.*cos(ups));
xb0   = r0.*cos(ups);
yb0   = r0.*sin(ups);

% negatives kappa3
ka_31 = -0.02*mu_G;
B1    = -3*ka_31*mu_G/L^4;
r1    = p./(1+ex.*cos(ups-B1*ups))*1./(1+(B1.*C*p./...
        (1+ex.*cos(ups-B1*ups))));
xb1   = r1.*cos(ups);
yb1   = r1.*sin(ups);

ka_32 = 0.02*mu_G;
% positives kappa3
B2    = -3*ka_32*mu_G/L^4;
r2    = p./(1+ex.*cos(ups-B2*ups))*1./(1+(B2.*C*p./...
        (1+ex.*cos(ups-B2*ups))));
xb2   = r2.*cos(ups);
yb2   = r2.*sin(ups);


%--------------------------------------------------------------------------
%%
% Graphische Ausgabe

header='Simulation von gestoerten Keplerbahnen mit Stoerterm kappa3';
figure('Name',header);
plot(xb1, yb1, 'Color', Colors(2,:),'LineWidth',1);
hold on;
plot(xb2, yb2, 'Color', Colors(4,:),'LineWidth',1);
plot(xb0, yb0, 'Color', Colors(3,:),'LineWidth',1);
axis equal;
ylim([-1.2 1.2]);
xlim([-1.2 1.2]);
xlabel('x')
ylabel('y');
ttlstr=strjoin(['\rm \it e\rm = ',string(num2str(ex)),...
    '\rm \it a\rm = ',string(num2str(a))]);
title(ttlstr, 'FontSize', 14);
legend(strjoin(['\kappa_3 = ',string(num2str(ka_31,'%+5.2f'))]),...
      strjoin(['\kappa_3 = ',string(num2str(ka_32,'%+ 5.2f'))]),...
      '\kappa_3 = 0');
legend box off;
grid on;
ax=gca;
set(gca,'Fontsize', 14, 'linewidth', 1);

header='Simulation von gestoerten Keplerbahnen mit Stoerterm kappa3';
figure('Name',header);
plot(ups/pi, r1./r0, 'Color', Colors(2,:),'LineWidth',1);
hold on;
plot(ups/pi, r2./r0, 'Color', Colors(4,:),'LineWidth',1);
plot(ups/pi, r0./r0, 'Color', Colors(3,:),'LineWidth',1);
ylim([0.85 1.15]);
xlim([0 4]);
xlabel('Wahre Anomalie {\upsilon} in \pi')
ylabel('\it r/r\rm_{Kepler}');
ttlstr=strjoin(['\rm \it e\rm = ',string(num2str(ex)),...
    '\rm \it a\rm = ',string(num2str(a))]);
title(ttlstr, 'FontSize', 14);
legend(strjoin(['\kappa_3 = ',string(num2str(ka_31,'%+5.2f'))]),...
      strjoin(['\kappa_3 = ',string(num2str(ka_32,'%+ 5.2f'))]),...
      '\kappa_3 = 0');
legend box off
hold on;
grid on;
ax=gca;
set(gca,'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------