% -------------------------------------------------------------------------
% KarussellCoriolis01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Karussell und Corioliskraft 
% Programm berechnet die Effekte der Corioliskraft auf einem Karussell.
% Die Reibung wird berücksichtigt.
% Numerische Lösung mit ODE45.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":","-.",":"];

%%
% Anfangsbedingungen, Parameter

T0     = 180;                       % Rotationsperiode in s
omega0 = 2*pi/T0;                   
omeg = omega0;
omegdot = 0;
maxzz = 5000;                       % Stuetzstellen
RK    = 50;                         % Radius Karussell
alpha = linspace(0,360);
KreisX = RK*cosd(alpha);
KreisY = RK*sind(alpha);

%%
% Numerisch ohne Reibung und ohne Beschleunigung
% Konstante Rotationsperiode T =180 S

% Fahrzeit ohne Reibung 3
v0    = 5;                          % Geschwindigkeit des Go-Cart in m/s
tmax = 3*2*RK/v0;
tv   = linspace(0,tmax,1000);
eta  = 0.0000;  % etwa 0
iend = 7;

% Konstante Geschwindigkeit, variabler Winkel
for k = 1:iend
    phi = (k-1)*15+90;                % Abfahrtswinkel in °
    phistr(k,:) = strcat("\phi=",string(num2str(phi,'%2d °')));
    v0x   = v0*cosd(phi);             % Geschwindigkeit x
    v0y   = v0*sind(phi);             % Geschwindigkeit y
    start = [0;v0x;-RK;v0y];
    options = odeset('RelTol',1e-6);     
    [t,Y]   = ode45(@CoriolisDGLS, tv, start, options, omeg, omegdot, eta);
    % Koordinaten
    Koord(k).t(:)    = t(:);
    Koord(k).xw(:)   = Y(:,1);
    Koord(k).yw(:)   = Y(:,3);
    Koord(k).rw(:)   = sqrt(Y(:,1).^2+Y(:,3).^2);
end
phistr(iend+1,:) = "Karussellrand";
Koord = CheckRand(Koord,iend,RK);

figure();
subplot(1,2,1);
for k=1:iend
   plot(Koord(k).xw(:),Koord(k).yw(:),'LineWidth',1,'Color',Colors(k,:),...
        'LineStyle',Style(k));
   hold on;
end
grid on;
xlabel('x in m ','FontSize',14);
ylabel('y in m ','FontSize',14);
axis square;
ylim([-RK,RK]);
xlim([-RK,RK]);
plot(KreisX,KreisY,'LineWidth',1,'Color',Colors(9,:),...
        'LineStyle',Style(1));
legend(phistr, 'location','southeast','NumColumns',2);
legend box off;
titlestr1= strcat("\it{v_0} \rm{= }",num2str(v0,'%3.1f m/s    '));
titlestr2= strcat("   {\eta} \rm{= }" ,num2str(eta,'%3.0f '));
titlestr3= strcat("   {T_0} \rm{= }" ,num2str(T0,'%3.0f s'));
title(titlestr1+titlestr2+titlestr3,'FontWeight','normal','FontSize',14);
set(gca,'Fontsize', 16);


% Konstanter Winkel, variable Geschwindigkeit 
for k = 1:iend
    phi = 135;
    v0  = k;
    phistr(k,:) = strcat("\it{v_0} =",string(num2str(v0,'%2d m/s')));
    v0x   = v0*cosd(phi);             % Geschwindigkeit x
    v0y   = v0*sind(phi);             % Geschwindigkeit y
    start = [0;v0x;-RK;v0y];
    options = odeset('RelTol',1e-6);     
    [t,Y]   = ode45(@CoriolisDGLS, tv, start, options, omeg, omegdot, eta);
    % Koordinaten
    Koord(k).t(:)    = t(:);
    Koord(k).xw(:)   = Y(:,1);
    Koord(k).yw(:)   = Y(:,3);
    Koord(k).rw(:)   = sqrt(Y(:,1).^2+Y(:,3).^2);
end
Koord = CheckRand(Koord,iend,RK);
subplot(1,2,2);
for k=1:iend
   plot(Koord(k).xw(:),Koord(k).yw(:),'LineWidth',1,'Color',Colors(k,:),...
        'LineStyle',Style(k));
   hold on;
end
grid on;
xlabel('x in m ','FontSize',14);
ylabel('y in m ','FontSize',14);
axis square;
ylim([-RK,RK]);
xlim([-RK,RK]);
plot(KreisX,KreisY,'LineWidth',1,'Color',Colors(9,:),...
        'LineStyle',Style(1));
legend(phistr, 'location','southeast','NumColumns',2);
legend box off;
titlestr1= strcat("\phi \rm{= }",num2str(phi,'%3.0f°   '));
titlestr2= strcat("   {\eta} \rm{= }" ,num2str(eta,'%3.0f  '));
titlestr3= strcat("   {T_0} \rm{= }" ,num2str(T0,'%3.0f s '));
title(titlestr1+titlestr2+titlestr3,'FontWeight','normal','FontSize',14);
set(gca,'Fontsize', 16);
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen

function dY = CoriolisDGLS(t, Y, Omeg, OmegDot, Eta)
    dY    = zeros(4,1);
    dY(1) = Y(2);
    dY(2) = Omeg^2*Y(1) + 2*Omeg*Y(4) + OmegDot*Y(3)- Eta*Y(2);
    dY(3) = Y(4);
    dY(4) = Omeg^2*Y(3) - 2*Omeg*Y(2) - OmegDot*Y(1)- Eta*Y(4);
end


function Koord = CheckRand(KoordAlt,IEnd, RK)
    Koord = KoordAlt;
    for iPlot = 1:IEnd
       Fac = 1;
       for k=2:length(KoordAlt(iPlot).rw(:))
       if abs(KoordAlt(iPlot).rw(k)) > RK*Fac 
          Fac = 0; 
          Koord(iPlot).xw(k) = NaN; 
          Koord(iPlot).yw(k) = NaN; 
          Koord(iPlot).rw(k) = NaN; 
       end
      end  
    end
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
