% -------------------------------------------------------------------------
% Gezeitenkraefte.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Wirkung von Gezeitenkräften auf einen Ring von 
% 36 Massenpunkten, die im Feld einer großen Masse frei fallen.
% 
%--------------------------------------------------------------------------


clear all
close all 
addpath('./IncludeFolder/')
Colors=GetColorLines;

% Vorbereitung Felder
alpha    = zeros(1,36);
phi      = zeros(1,36);
R        = 2;
r0       = 10*R;
M        = 1;
G        = 1;

k   = 1:36;
phi = k*10;
salpha = R*sind(phi)./sqrt(1-2*R*cosd(phi)/r0+(R/r0)^2)/r0;
calpha = sqrt(1-salpha.^2);

xi0 = R*sind(phi);
zi0 = r0-R*cosd(phi);

TSPAN = linspace(0,75,101);
options = odeset('AbsTol',1.e-9,'RelTol',1.e-7);

% MATLAB Runge-Kutta (4,5) Ode45 Verfahren
for k=1:length(xi0)
    AB    =[xi0(k),0,zi0(k),0];           % Anfangsbedingungen
    [tout,Y] = ode45(@Diff_EquN,TSPAN,AB,[],salpha(k),calpha(k),M*G);  
%   t(:,k)     = tout(:);
    x(:,k)     = Y(:,1);   
    z(:,k)     = Y(:,3);   
end

AB    =[0,0,r0,0];                      % Anfangsbedingungen
[tout,Y] = ode45(@Diff_EquN,TSPAN,AB,[],salpha(k),calpha(k),M*G);  
% t(:,37)     = tout(:);
x(:,37)     = Y(:,1);   
z(:,37)     = Y(:,3);   

% Ausgabe Linien konstanter Deklination
figure('Name',' Verformung Ring im Gravitationsfeld');
hold on
plot(xi0,zi0,'ro');
plot(0,r0,'ro');
plot(x(end,:),z(end,:),'gd');
PlotCircle (0,0,0.5,'k',2);             % große Masse
axis equal;
grid on;
axis equal
ylim([-5 25]);
xlim([-5 5]);
xlabel('a.u.');
ylabel('a.u.');
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Funktionen

% DGL System
% Differentialgleichungssystem
function  dY  =  Diff_EquN(t, Y, salpha, calpha, MG)
    dY  =  zeros(4,1);  % es muss ein Spaltenvektor zurückgegeben werden 
    dY(1)  =  Y(2);
    dY(2)  = -MG*salpha/(Y(1)^2+Y(3)^2);
    dY(3)  =  Y(4);
    dY(4)  = -MG*calpha/(Y(1)^2+Y(3)^2);
end


% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------

