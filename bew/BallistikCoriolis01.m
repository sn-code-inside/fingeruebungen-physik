% -------------------------------------------------------------------------
% BallistikCoriolis01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Ballistik und Corioliskraft 
% 
% Programm berechnet die Effekte der Corioliskraft auf den freien Fall.
% Exakte Lösung und bei Berücksichtigungd der Reibung numerisache Lösung 
% mit ODE45.
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
% Aufgabenteil 1
% Exakt und numerisch bei Reibung 
% Achtung x nach Süden 

% Wurfzeit ohne Reibung *1.1
hz   = 1000;    
tmax = 2*sqrt(2*hz/g)
tv   = linspace(0,tmax,100);
c_w  = 0.45;
rho_L= 1.225;       % Dichte von Luft kg/m^3
% rho_K= 7900;      % Dichte von Eisen kg/m^3
rho_K= 790.0;       % Dichte von Holz kg/m^3
R    = 0.1;
eta0 = 3*c_w*rho_L/rho_K/R/8;
eta  = eta0;
omeg = 1.0*omeg;
iend = 4;
% Anfangsbedingungen Freier Fall
start = [0;0;0;0;hz;0];

for k = 1:iend
    % exakte Lösung
     % Näherungslösung
end

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
    xe(k,:)  = (g*sind(phi)*cosd(phi)/4/omegaE0^2).*(2*(omegaE0*tv).^2-...
         (1-cos(2*omegaE0*tv)));
    ye(k,:)  = (g*cosd(phi)/2/omegaE0).*(tv-sin(2*omegaE0.*tv)/2/omegaE0);
    ze(k,:)  = hz - (g*tv.^2)/2 + g*cos(phi)*cosd(phi)/4/omegaE0/...
         omegaE0*(2*(omegaE0*tv).^2-(1-cos(2*omegaE0*tv)));

    Koord(k).t(:)    = t(:);
    Koord(k).zw(:)   = Y(:,5);
    Koord(k).xw(:)   = Y(:,1);
    Koord(k).yw(:)   = Y(:,3);
end

for iPlot = 1:iend*2
  for k=1:length(Koord(iPlot).zw(:))
   if Koord(iPlot).zw(k) < -25
      Koord(iPlot).zw(k) = NaN; 
      Koord(iPlot).xw(k) = NaN; 
      Koord(iPlot).yw(k) = NaN; 
   end
  end  
end

for iPlot = 1:iend*2
  for k=1:length(ze(iPlot,:))
   if ze(iPlot,k) < -10
      ze(iPlot,k) = NaN; 
      xe(iPlot,k) = NaN;
      ye(iPlot,k) = NaN;
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
   plot(Koord(k).yw(:),Koord(k).zw(:),'LineWidth',1,'Color',Colors(cp,:),...
        'LineStyle',Style(kz));
   hold on;
end
grid on;
xlabel('Ostabweichung in m ','FontSize',14);
ylabel('Höhe in m ','FontSize',14);
ylim([0,inf]);
xlim([0,inf]);
legend(phistr, 'location','north','NumColumns',2);
legend box off;
titlestr1= strcat("\it{h} \rm{=}",num2str(hz,' %4d m    '));
titlestr2= strcat(" {\eta} \rm{=}" ,num2str(eta0,' %3.2e 1/m'));
title(titlestr1+titlestr2,'FontWeight','normal','FontSize',14);
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
   plot(Koord(k).t(:),Koord(k).yw(:),'LineWidth',1,'Color',Colors(cp,:),...
       'LineStyle',Style(kz));
   hold on;
end
% Zum Vergleich exakte Lösung
% for k= 5:iend*2
%   plot(tv(:),ye(k,:),'LineWidth',2,'Color',Colors(10,:),...
%        'LineStyle','--');
%   hold on;
% end
grid on;
ylabel('Ostabweichung in m ','FontSize',14);
xlabel('Fallzeit in s ','FontSize',14);
xlim([0,inf]);
legend(phistr, 'location','northwest','NumColumns',2);
legend box off;
titlestr1= strcat("\it{h} \rm{=}",num2str(hz,' %4d m    '));
titlestr2= strcat(" {\eta} \rm{=}" ,num2str(eta0,' %3.2e 1/m'));
title(titlestr1+titlestr2,'FontWeight','normal','FontSize',14);

set(gca,'Fontsize', 16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen
% Differentialgleichungssystem für Corioliskraft
function dY = CoriolisDGLS(t,Y,Omeg,Cphi,Sphi,G,Eta)
    vbet  = sqrt(Y(2)^2+Y(4)^2+Y(6)^2);
    dY    = zeros(6,1);
    dY(1) = Y(2);
    dY(2) = 2*Omeg*Sphi*Y(4)-Eta*Y(2)*vbet;
    dY(3) = Y(4);
    dY(4) = -2*Omeg*Sphi*Y(2)-2*Omeg*Cphi*Y(6)-Eta*Y(4)*vbet;
    dY(5) = Y(6);
    dY(6) = 2*Omeg*Cphi*Y(4)-G-Eta*Y(6)*vbet;  
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
