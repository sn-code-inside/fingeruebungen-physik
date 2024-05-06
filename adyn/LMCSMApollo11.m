%% -------------------------------------------------------------------------
% LMCSMApollo11.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Rendezvous von LM und CSM von Apollo 11 in der Mondumlaufbahn.
% Berechnung der Kinematik des Weltall-Rendezvous über die CW-Gleichungen.
% -------------------------------------------------------------------------

%%
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-.", "-", ":", "--", ":"];

%% Parameter Numerische Lösung
% 
h  = 111.12;            % Höhe des CSM über Mondoberfläche in km
Delh = 27;              % Differenzhöhe von LEM und CSM in km
MM = 7.3483e22;         % Masse Mond in kg
RM = 1.738e6;           % Mondradius in m
G  = 6.671e-11;         % G in (m^3 / kg /s^2)
RCSM = RM+h*1e3;        % Umlaufradius des CSM (Target) in m
RLEM = RCSM-Delh*1e3;   % Umlaufradius desLEM in m

wCSM = sqrt(G*MM/RCSM^3); % Umlaufperiode des CSM (Target) in s
wLEM = sqrt(G*MM/RLEM^3); % Umlaufperiode des LM (Jäger) in s

TCSM = 2*pi/wCSM;    
T0min= TCSM/60;         % Umlaufperiode des CSM (Target) in min

tmax  = 42*60;          % Simulationszeit in s
tspan = linspace(0,tmax,42+1);
t= tspan;

%% Numerische Auswertung CW-Gleichung 

% Anfangswerte (Position und Geschwinddigkeit ISS(Target)) in m,m/s
x= wrapTo360(480)

w0      = wCSM;
tflight = tmax;         % vorgegebene Flugzeit von TPI zu TPF
for k= 1:3
    y0(k)   = -27.78e3+(k-2)*2000; % y-Anfansgposition LEM in m im CSM-KOS
    x0      = y0(k)*cotd(26.5);    % x-Anfansgposition LEM in m im CSM-KOS
    vLEM0x  = -1.5*w0*y0(k);       % Geschwindigkeit des LEM vor TPI Burn
    vLEM0y  = 0;

    % Bestimmung Anfangsgeschwindigkeiten
    sinw0 = sin(w0*tflight); cosw0 = cos(w0*tflight);
    C       = 3*w0*tflight*sinw0-8*(1-cosw0);
    dx0      = (w0*x0*sinw0-w0*y0(k)*(6*w0*tflight*sinw0-14*(1-cosw0)))/C;
    dy0      = (-2*w0*x0*(1-cosw0)+w0*y0(k)*(4*sinw0-3*w0*tflight*cosw0))/C;
    v0       = sqrt(dx0^2+dy0^2);
    gamma    = atand((dy0-vLEM0y)/(dx0-vLEM0x));
    lgdstr2(k,:)=num2str(y0(k)/1e3,' y_0 = %3.1f km  ');
    lgdstr3(k,:)=num2str(sqrt((dy0-vLEM0y)^2+(dx0-vLEM0x)^2),...
                 ' Delta v = %3.1f m/s');
    lgdstr4(k,:)=num2str(gamma,    ' gamma = %3.1f ° ');

    xCW(k,:)  = (6*y0(k)+4*dx0/w0)*sin(w0*t)+2*dy0*cos(w0*t)/w0-2*dy0/w0 ...
            -(3*dx0+6*w0*y0(k))*t + x0;
    yCW(k,:)  = -(3*y0(k)+2*dx0/w0)*cos(w0*t)+dy0*sin(w0*t)/w0...
             +4*y0(k)+2*dx0/w0;
end

%% Graphische Ausgabe

figure(1)
hold on
for k=1:3
    hp(k) = plot(xCW(k,:)/1e3, yCW(k,:)/1e3,'+','MarkerIndices',...
            1:2:length(xCW),'color',Colors(k,:),'MarkerSize',6,...
            'LineWidth',1,'LineStyle',Style(k));
    plot(xCW(k,1)/1e3,y0(k)/1e3,'d','MarkerSize',10,'color',...
            Colors(k,:),'LineWidth',2);

end
plot(0,0,'s','MarkerSize',10,'color',Colors(4,:),'LineWidth',2);
xL = [xCW(1,1) 0]/1e3;
yL = [y0(1) 0]/1e3;
line(xL,yL,'Color','red','LineStyle','--')
axis equal
axis ([x0 10e3 y0(1) 0]*1.2e-3)          %windows size
ylabel ('\it y \rm in km')
xlabel ('\it x \rm in km')
grid on 
legend(hp,lgdstr3,'location','southeast','Numcolumns',1)
legend box off
h=title('Kopplung Apollo 11 LEM und CSM');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
