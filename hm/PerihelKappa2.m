% -------------------------------------------------------------------------
% PerihelKappa2.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Perihelverschiebung mit Störterm kappa2 ungleich Null (vgl. Kap. 2.5.4).
% Die Parameter wurden höher gewählt als sie in der Realität
% sind, um die Auswirkungen graphisch besser zu veranschaulichen. 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
ups0 = 0;
ups=linspace(0,8*pi, 4000); % vier Perioden

%Parameter
a=1;                        % Radius
mu_G=1;                     % Gravitationskonstante * Masse 
ex=0.5;                     % Exzentrizitaet 
ka_2 = -0.2*mu_G;           % kappa2 

ka_2v(1) = ka_2/2;          % Erzeugen verschiedener Werte von kappa2 
ka_2v(2) = -ka_2/2;
ka_2v(3) =  ka_2;
ka_2v(4) = -ka_2;

for k =1:4
    % Erzeugung der Ueberschriften 
    ttlstr(k,:)=strjoin(['\rm \it e\rm = ',string(num2str(ex)),...
    '\rm \it a\rm = ',string(num2str(a)),' {\kappa}_2 = ', ...
           string(num2str(ka_2v(k)))]);
end

%--------------------------------------------------------------------------
%%
%Graphische Ausgabe
header='Simulation von gestoerten Keplerbahnen mit Stoerterm kappa2';
figure('Name',header);

% Iteration durch verschiedene Werte von kappa2
for k=1:4
    subplot(2,2,k)

    ka_2 = ka_2v(k);
    % Berechnung des Drehimpuls 
    L=sqrt(a*mu_G*(1-ex*ex));
    % Berechnung der Energie (auf die Masse normiert)
    H= -(mu_G*mu_G*(1-ex*ex))/2/L/L;

    % Berechnung der Parameter aus Formel 3.131 
    pk  = L*L*(1+2*ka_2/L/L)/mu_G;
    exk = ex*sqrt(1+2*H*ka_2/(mu_G*mu_G+2*H*L*L));
    Ck  = -(1+2*ka_2/L/L);

    r   = pk./(1+exk.*cos(Ck*(ups-ups0)));
    % Berechnung fuer die Darstellung im kartesischen KOS
    xb   = r.*cos(ups);
    yb   = r.*sin(ups);

    plot(xb, yb, 'Color', Colors(k,:),'LineWidth',1);
    ylim([-2 2]);
    xlim([-2 2]);
    hold on;
    grid on;
    grid minor;
    xlabel('x')
    ylabel('y');
    title(ttlstr(k,:), 'FontSize', 14);
    axis equal;
    ax=gca;
    set(gca,'Fontsize', 14, 'linewidth', 1);
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
