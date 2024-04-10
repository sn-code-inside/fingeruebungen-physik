% -------------------------------------------------------------------------
% Expander.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Dynamik des Nichtlinearen Pendels
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Symbolische Berechenung

% Definiert die verallgemeinerten Koordinaten theta1, theta2 und Parameter

syms m k L0 qmax x dx g eta 'real'

% Verallgemeinerte Koordinate und Ableitung
q  = x;
dq = dx;

% kinetische Energie
T_trans = 1/2* m *(dx*dx);
T_rot = 0;
T_kin = T_trans + T_rot;

% potentielle Energien
U_pot =  k* L0^2*(sqrt(1+(x/L0)^2)-1)^2;

% Lagrange-Funktion
L = T_kin - U_pot;

%% Ausgabe der Ergebnisse
fprintf("------------------------------------\nLagrange-Funktion L:");
fprintf("\n------------------------------------\n");
fprintf(1, 'T = %s \n', T_kin);
fprintf(1, 'U = %s \n', U_pot);
fprintf(1, 'L = T - U = %s \n', L);
fprintf("\n");

% Externe Kräfte 
Q = 0;

% Berechnet Euler-Lagrange-Gleichungen:
EQ = EulerLagrange(q,dq,L,Q,true);


%%
NPts   = 500;
m      = 2;
k      = 1000;
L0     = 2;
tmax   = 10;           % Zeitspanne
tspan  = linspace(0,tmax,NPts);
dx0    = 0;

% Numerische Lösung
for iplot = 1:4
    qmax    = 0.05+0.05*iplot;         % Anfangsamplitude 
    AB=[qmax;0];       % AB für ode45
    opt=odeset('AbsTol',1.e-8,'RelTol',1.e-5);      
    [t,Y]=ode45(@dgl_expander,tspan,AB,opt,k,L0,m); % Numerische Lösung
    xN(:,iplot) = Y(:,1);
    
    % Analytische Näherung Oszillator mit modifizierter Frequenz
    T(iplot)     = 2*sqrt(2*pi)*sqrt(m*4*L0^2/k)*gamma(1.25)/...
                   gamma(0.75)/qmax;
    omega(iplot) = 2*pi/T(iplot);
    xH(:,iplot) = qmax*cos(omega(iplot)*tspan(:));
    lgdstr(iplot,:) = num2str(qmax/L0,'%5.3f');
end
%% Graphische Darstellung

figure()
hold on;
for iplot = 1:4
   h(iplot)= plot(t,xN(:,iplot)/L0,'Color',Colors(iplot,:),'LineWidth',1,...
    'LineStyle', Style(1));
   plot(tspan,xH(:,iplot)/L0,'Color',Colors(iplot,:),'LineWidth',2,...
    'LineStyle', Style(3));
end
hl=legend(h,lgdstr,'NumColumns',4);
set(hl,'FontSize',12)

axis([0 max(tspan) -qmax*1.25/L0 qmax*1.25/L0 ]);
xlabel('Zeit in s','FontSize',14)
ylabel('Amplitude q/L_0 ','FontSize',14)
grid on;
legend box off;
set(gca,'FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%

function dY = dgl_expander(t, Y, k, L0, m)
    x  = Y(1);
    dY = [Y(2); (2*k*x)/(m*(x^2/L0^2 + 1)^(1/2)) - (2*k*x)/m];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
