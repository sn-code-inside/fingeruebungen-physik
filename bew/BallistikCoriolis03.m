% -------------------------------------------------------------------------
% BallistikCoriolis01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Effekte der Corioliskraft auf die Ballistik
% Hier Schuss eines Projektils nach Süden. Die Reibung wird berücksichtigt.
% Numerische Lösung mit ODE45.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Anfangsbedingungen, Parameter
g = 9.81;                           % Schwerebeschleunigung
% Kreisrotationsfrequenz Erde
omegaE0 = 2*pi/86164;  
maxz = 5000;
omeg = omegaE0;
v0    = 400;                        % Geschwindigkeit in m/s

%%
% Aufgabenteil 3
% Numerisch mit Reibung 
% Achtung x nach Süden 

AZ    =  0;                         % Azimut oder Vorhaltewinkel
alpha = 10;                         % Abschußwinkel in °
v0x   = v0*cosd(alpha)*cosd(AZ);    % Schußrichtung Süden
v0y   = v0*cosd(alpha)*sind(AZ);    % Richtung Osten
v0z   = v0*sind(alpha);

% Wurfzeit ohne Reibung *1.1
tmax = 1.5*2*v0z/g;
tv   = linspace(0,tmax,1000);
c_w  = 0.45;                        % Widerstandsbeiwert
rho_L= 1.225;                       % Dichte Luft in kg/m^3
rho_K= 7900;                        % Dichte Eisen in kg/m^3
R    = 0.1;                         % Radius Kanonenkugel
eta0 = 3*c_w*rho_L/rho_K/R/8;
eta  = eta0;
omeg = 1.0*omeg;
iend = 4;
hz   = 500;    
% Anfangsbedingungen Ballistik
start = [0;v0x;0;v0y;0;v0z];
for k = 1:iend*2
    kz = k;
    if k > iend 
        kz = k - iend;
        eta = 0;
    end
    phi = (kz-1)*20+10;
    if k<5 phistr(k,:) = strcat("\phi=",string(num2str(phi,'%2d °')));
    else phistr(k,:) = "o. Reibung"; end

    cphi = cosd(phi);
    sphi = sind(phi);
    options = odeset('RelTol',1e-6);     
    [t,Y]   = ode45(@CoriolisDGLS, tv, start, options, omeg, cphi,...
                    sphi, g, eta);
    % Koordinaten
    Koord(k).t(:)    = t(:);
    Koord(k).zw(:)   = Y(:,5);
    Koord(k).xw(:)   = Y(:,1);
    Koord(k).yw(:)   = Y(:,3);
end

for iPlot = 1:iend*2
  for k=1:length(Koord(iPlot).zw(:))
   if Koord(iPlot).zw(k) < -5
      Koord(iPlot).zw(k) = NaN; 
      Koord(iPlot).xw(k) = NaN; 
      Koord(iPlot).yw(k) = NaN; 
   end
  end  
end

figure();
subplot(1,2,1);
for k=1:iend*2
   if k > iend 
       kz = 5;
       cp = k -iend;
   else
       kz = 1;
       cp = k;
   end
   plot(Koord(k).xw(:),Koord(k).zw(:),'LineWidth',1,'Color',Colors(cp,:),...
        'LineStyle',Style(kz));
   hold on;
end
grid on;
xlabel('Weite in m ','FontSize',14);
ylabel('Höhe in m ','FontSize',14);
ylim([0,500]);
xlim([0,inf]);
legend(phistr, 'location','north','NumColumns',2);
legend box off;
titlestr1= strcat("\it{v_0} \rm{=}",num2str(v0,' %3.0f m/s    '));
titlestr2= strcat(" {\eta} \rm{=}" ,num2str(eta0,' %3.2e 1/m'));
titlestr3= strcat(" {\alpha} \rm{=}" ,num2str(alpha,' %3.1f°'));
title(titlestr1+titlestr2+titlestr3,'FontWeight','normal','FontSize',14);

set(gca,'Fontsize', 16);

subplot(1,2,2);
for k=1:iend*2
   if k > iend 
       kz = 5;
       cp = k -iend;
   else
       kz = 1;
       cp = k;
   end
   plot(Koord(k).xw(:),Koord(k).yw(:),'LineWidth',1,'Color',Colors(cp,:),...
       'LineStyle',Style(kz));
   hold on;
end
grid on;
ylabel('Ostabweichung in m ','FontSize',14);
xlabel('Weite in m ','FontSize',14);
xlim([0,inf]);
legend(phistr, 'location','south','NumColumns',2);
legend box off;
title(titlestr1+titlestr2+titlestr3,'FontWeight','normal','FontSize',14);
set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen

function dY = CoriolisDGLS(t,Y, Omeg, Cphi, Sphi, G, Eta)
    vbet  = sqrt(Y(2)^2+Y(4)^2+Y(6)^2)
    dY    = zeros(6,1);
    dY(1) = Y(2);
    dY(2) = 2*Omeg*Sphi*Y(4) - Eta*Y(2)*vbet;
    dY(3) = Y(4);
    dY(4) = -2*Omeg*Sphi*Y(2) - 2*Omeg*Cphi*Y(6) - Eta*Y(4)*vbet;
    dY(5) = Y(6);
    dY(6) = 2*Omeg*Cphi*Y(4) - G - Eta*Y(6)*vbet;  
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

