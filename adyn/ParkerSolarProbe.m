% -------------------------------------------------------------------------
% ParkerSolarProbe.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet verschieden Gravity-Assist-Manöver der Parker-Solar-
% Probe und simuliert die Bahn im Zeitraum von Start bis 2024.
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style  = ["-", "-.", ":", "--", ":"];
NamesP = ["Merkur ","Venus  ","Erde   ","Mars   ","Jupiter",...
             "Saturn ","Uranus ","Neptun"];

% Parameter, hier alles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % m
aV  = 0.72333199*AE;  % Bahnradius Venus in m
aV  = 0.723*AE;       % Bahnradius Venus in m
MS  = 1.989e30;       % Masse Sonne in kg 
MV  = 4.867e24;       % Masse Venus in kg

% Planetendaten, ab hier alles in km,s,kg
GS   = MS*G*1e-9;   % Sonne
GV   = MV*G*1e-9;   % Venus
RV   = 6051;        % Planetenradius Venus in km
R0   = aV/1000;     % Radius SOI in km
v0   = sqrt(GS/R0);         % Umlaufgeschwindigkeit Venus @ SOI

% EVENTS
%   2018
%     2018-Oct-03: Venus flyby #1 @ 8:44 UTC
%   2019
%     2019-Dec-26: Venus flyby #2
%   2020
%     2020-Jul-11: Venus flyby #3
%   2021
%     2021-Feb-20: Venus flyby #4
%     2021-Oct-16: Venus flyby #5
%   2023
%     2023-Aug-21: Venus flyby #6
%   2024
%     2024-Nov-06: Venus flyby #7 (last one)
%     2024-Dec-24: Perihelion #22 (first close approach)
%     
    
% Perihel in km
rP   = [0.207 0.166 0.130 0.095 0.074 0.062 0.053, 0.046]*AE/1000; 
% Aphel in km
rA   = [1.013 0.938 0.874 0.817 0.783 0.761 0.745 0.731]*AE/1000; 
a    = 0.5*(rA+rP);
ecc  = 1-rP./a;
% Geschwindigkeit Sonne im Kreuzungspunkt
vR   = v0*sqrt(2-R0./a);
vt   = v0*sqrt(a.*(1-ecc.^2)/R0);
vr   = sqrt(vR.^2-vt.^2);
vinf = sqrt(vr.^2+(vt-v0).^2);
% Umlaufzeit
Tp   = 2*pi./sqrt(GS./a.^3);   %in s
Tpd  = Tp/86400;  %in Tagen
% Perihelgeschwindigkeit
vP   = sqrt(GS*(2./rP-1./a));
% Winkel zwischen vinf und vP
beta1d = acosd((v0-vt)./vinf);

%% Berechnung maximal möglicher Streuwinkel und Winkel zw. vinf und vP
%  (alles in km /s)
vF   = sqrt(2*GV./(RV+500)); % Perizentrums-Fluchtgeschwindigkeit in km/s
% Winkel zwischen vinf und vP
thetad = 2*asind(1./(1+2*vinf.^2./vF.^2));
totalthetad = sum(thetad);

%% Berechnung realer Streuwinkel und Winkel zw. vinf und vP
%  (alles in km /s)
RVreal = [1.67 1.49 1.15 1.43 1.70 1.79 1.37 100000]*RV;
vFreal   = sqrt(2*GV./RVreal); % Perizentrums-Fluchtgeschwindigkeit in km/s
% Winkel zwischen vinf und vP
thetad_real = 2*asind(1./(1+2*vinf.^2./vFreal.^2));
totalthetad_real = sum(thetad_real);

%% Ausgabe Werte

nrorbits = ["a" "b" "c" "d" "e" "f" "g" "h"];
fprintf('\n Orbit |  r_P  (AE)  | r_A  (AE)  |  beta_1 (°) |    e     |   T_p (Tage) | v_p (km/s)\n');
for k=1:length(vP)
    fprintf('  %s    |   %5.3f     |   %5.3f    |   %5.2f     |  %5.3f   |   %6.1f     |  %6.2f \n',...
            nrorbits(k), rP(k)*1000/AE, rA(k)*1000/AE, beta1d(k), ecc(k), Tpd(k), vP(k));
end
fprintf('\n');
fprintf('\n GA-Manöver  |  theta (°)  |  theta_total (°)  | R_V,real/R_V ');
for k=1:length(vP)-1
    fprintf('\n %u  %s -> %s   |    %5.3f    |     %6.3f        |  %5.2f  ',...
            k, nrorbits(k), nrorbits(k+1), thetad_real(k), sum(thetad_real(1:k)), RVreal(k)/RV);
end
fprintf('\n');

%% Graphische Ausgabe

% Laden der Ephemeriden der Planeten und der PSP
Venus = importfile('Venus.csv');
VenusFlyBy= importfile('VenusFlyBy.csv'); %FlyBy Daten
Erde  = importfile('Erde.csv');
PSP   = importfile('PSP.csv');
xE  = Erde.x;   yE= Erde.y;
xV  = Venus.x;  yV= Venus.y;
xVFlyBy = VenusFlyBy.x;  yVFlyBy = VenusFlyBy.y;
xP = PSP.x;    yP= PSP.y;  %PSP

%Darstellung Kreisbahn Erde, Venus (Näherung)
u   = linspace(0,360,361);
xKE = Erde.r(1)*cosd(u); yKE = Erde.r(1)*sind(u);
xKV = aV*cosd(u)/AE; yKV = aV*sind(u)/AE;

% Trajektorie
figure()
hold on
axis equal;
axis square;
grid on
PlotCircle(0,0,0.02,Colors(10,:),2);
plot(xKV,yKV,'color',Colors(2,:),'LineWidth',1);
plot(xKE,yKE,'color',Colors(3,:),'LineWidth',1);
plot(xVFlyBy(1:7),yVFlyBy(1:7),'d','color',Colors(2,:),'markersize',6,...
    'linewidth',1,'markerfacecolor',Colors(2,:));
plot(xE(1),yE(1),'o','color',Colors(3,:),'markersize',7,...
    'markersize',10,'linewidth',1,'markerfacecolor',Colors(3,:));
plot(xP(1),yP(1),'d','color',Colors(4,:),'markersize',6,...
    'linewidth',1,'markerfacecolor',Colors(4,:));
plot(xP,yP,'color',Colors(4,:));
xlim([-0.3 1.1]);
ylim([-1.1 0.3]);
xlabel('x (AE)'); ylabel('y (AE)');
for kF=1:7
    if mod(kF,2) == 0
      text(xVFlyBy(kF)-0.05,yVFlyBy(kF),num2str(kF));
    else
      text(xVFlyBy(kF)+0.025,yVFlyBy(kF),num2str(kF));
    end
end
set(gca,'FontSize',14);
ttl=title("Parker-Solar-Probe");
set(ttl,'FontSize',14,'FontWeight','normal');


%Simulation

figure()
kstart =1;
for kF = 1:7 %length(VenusFlyBy.r)
    clf;
    hold on
    axis equal;
    axis square;
    grid on
    xlim([-0.3 1.1]);
    ylim([-1.1 0.3]);
    xlim([-1.1 1.1]);
    ylim([-1.1 1.1]);
    xlabel('x (AE)'); ylabel('y (AE)');
    set(gca,'FontSize',14);
    ttl=title("Parker-Solar-Probe");
    set(ttl,'FontSize',14,'FontWeight','normal');
    PlotCircle(0,0,0.02,Colors(10,:),2);
    plot(xKV,yKV,'color',Colors(2,:),'LineWidth',1);
    plot(xKE,yKE,'color',Colors(3,:),'LineWidth',1);
    if kF>1 
        plot(xVFlyBy(1:kF-1),yVFlyBy(1:kF-1),'d',...
            'color',Colors(2,:),'markersize',8,...
       'linewidth',1,'markerfacecolor',Colors(2,:)); 
    end
    if kF <2 
        plot(xE(1),yE(1),'o','color',Colors(3,:),'markersize',7,...
        'markersize',10,'linewidth',1,'markerfacecolor',Colors(3,:));
        plot(xP(1),yP(1),'d','color',Colors(4,:),'markersize',6,...
        'linewidth',1,'markerfacecolor',Colors(4,:));
    end
    while datetime(PSP.Datum(k)) < datetime(VenusFlyBy.Datum(kF)) 
        plot(xP(kstart:k),yP(kstart:k),'color',Colors(4,:));
        hP = plot(xP(k),yP(k),'s','color',Colors(4,:));
        hV = plot(xV(k),yV(k),'d','color',Colors(2,:),...
              'markersize',6,'linewidth',2);
        hT = text(0.5,1,strcat('t = ', datestr(datetime(PSP.Datum(k)))));
        pause(0.02);
        k = k +1;
        hV.Visible = 'off';
        hP.Visible = 'off';
        hT.Visible = 'off';
    end
    kstart = k;
    text(0.5,1,strcat('t = ', datestr(datetime(PSP.Datum(k)))));
    text(0.5,-1,strcat('FlyBy Nr:',num2str(kF)),...
        'FontSize',14,'FontWeight','normal','Color',Colors(2,:));
    plot(xVFlyBy(1:kF),yVFlyBy(1:kF),'d',...
            'color',Colors(2,:),'markersize',8,...
       'linewidth',1,'markerfacecolor',Colors(2,:)); 
    pause(2.5);
end
clf
hold on
axis equal;
axis square;
grid on
PlotCircle(0,0,0.02,Colors(10,:),2);
plot(xKV,yKV,'color',Colors(2,:),'LineWidth',1);
plot(xKE,yKE,'color',Colors(3,:),'LineWidth',1);
plot(xVFlyBy(1:7),yVFlyBy(1:7),'d','color',Colors(2,:),'markersize',6,...
    'linewidth',1,'markerfacecolor',Colors(2,:));
plot(xP,yP,'color',Colors(4,:));
xlim([-1.1 1.1]);
ylim([-1.1 1.1]);
xlabel('x (AE)'); ylabel('y (AE)');
for kF=1:7
    if mod(kF,2) == 0
      text(xVFlyBy(kF)-0.05,yVFlyBy(kF),num2str(kF));
    else
      text(xVFlyBy(kF)+0.025,yVFlyBy(kF),num2str(kF));
    end
end
set(gca,'FontSize',14);
ttl=title("Parker-Solar-Probe");
set(ttl,'FontSize',14,'FontWeight','normal');



%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

function Planet = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  Planet = IMPORTFILE(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  Planet = IMPORTFILE(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. 
%
%  Example:
%  Planet = importfile("...\Data\Venus.csv", [2, Inf]);
%


%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["Datum", "Laenge", "Breite", "Abstand", "Laenge1", "Breite1", "x", "y", "r"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "InputFormat", "yyyy-MM-dd");
opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Planet = readtable(filename, opts);

end
%Ende Funktionen
% -------------------------------------------------------------------------

