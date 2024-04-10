% -------------------------------------------------------------------------
% BallistikCoriolis02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Ballistik und Corioliskraft 
% 
% Programm berechnet die Effekte der Corioliskraft auf die Ballistik
% Exakte Lösung. Hier Schuss eines Projektils nach Süden.
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
for k=1:4
    phi(k) = (k-1)*20 + 10;
    phistr(k,:) = strcat("\phi=",string(num2str(phi(k),'%2d °')));
end

omeg = omegaE0;
v0    = 300;                        % Geschwindigkeit in m/s

%%
% Aufgabenteil 2 a
% Exakt ohne Reibung 
% Achtung x nach Süden 

AZ    = 0;                          % Abschußwinkel bzgl. x-Achse
alpha = 45;                         % Abschußwinkel in °
v0x   = v0*cosd(alpha)*cosd(AZ);    % Schußrichtung Süden
v0y   = v0*cosd(alpha)*sind(AZ);    % Richtung Osten
v0z   = v0*sind(alpha);
% Wurfzeit und Wurfweite ohne Reibung
tmax = 2*v0z/g*1.1;
t = linspace(0,tmax,maxz);
xmax = v0^2*sind(2*alpha)*cosd(AZ)/g;
% Wurfzeit ohne Reibung
t = linspace(0,tmax,maxz);

for k = 1:4
    phi = (k-1)*20+10;
    phistr(k,:) = strcat("\phi=",string(num2str(phi,'%2d °')));
    cphi = cosd(phi);
    sphi = sind(phi);
    A   = -(v0x*sphi+v0z*cphi);
    B   = v0y-g*cphi/2/omeg;
    arg = 2*omeg*t;
    % Höhe
    zw(k,:) = v0z*t -0.5*g*t.^2 + (A*(t-sin(arg)/2/omeg) + ...
              B*(1-cos(arg))/2/omeg)*cphi + 0.5*g*cphi^2*t.^2;
    % Ost-West-Richtung
    yw(k,:)  = v0y*t+(A*(1-cos(arg)) - B*(2*omeg*t-sin(arg)))/2/omeg; 
    % Süd-Nord-Richtung
    xw(k,:)   = v0x*t+(A*(t-sin(arg)/2/omeg)+B*(1-cos(arg))/2/omeg)*sphi+...
                0.5*g*sphi*cphi*t.^2;
end

phistr(5,:) = "Ohne Coriolis";
zw(5,:) = v0z*t -0.5*g*t.^2;
xw(5,:) = v0x*t; 
yw(5,:) = 0; 
 
for iPlot = 1:5
  for k=1:length(zw(iPlot,:))
   if zw(iPlot,k) < -5
      zw(iPlot,k)= NaN; 
      xw(iPlot,k)= NaN; 
      yw(iPlot,k)= NaN; 
   end
  end  
end


figure();
subplot(2,1,1);
for k=1:4 
    plot(xw(k,:),zw(k,:),'LineWidth',1,'Color',Colors(k,:),'LineStyle',Style(k));
    hold on;
end
plot(xw(5,:),zw(5,:),'LineWidth',1,'Color',Colors(9,:),'LineStyle','-');
grid on;
grid on;
xlabel('Weite in m ','FontSize',14);
ylabel('Höhe in m ','FontSize',14);
ylim([0,2500]);
xlim([0,9000]);
legend(phistr, 'location','south');
legend box off;
set(gca,'Fontsize', 16);

subplot(2,1,2);
for k=1:4 
    plot(xw(k,:),yw(k,:),'LineWidth',2,'Color',Colors(k,:),'LineStyle',Style(k));
    hold on;
end
plot(xw(5,:),yw(5,:),'LineWidth',2,'Color',Colors(9,:),'LineStyle',Style(5));
grid on;
grid on;
ylabel('Ostabweichung in m ','FontSize',14);
xlabel('Weite in m ','FontSize',14);
% ylim([0,inf]);
xlim([0,inf]);
legend(phistr, 'location','southwest');
legend box off;
set(gca,'Fontsize', 16);

%%
% Aufgabenteil 2 b
% Optimierung Abschuß auf Ziel in Schußweite auf 8000 m
% x nach Süden y nach Osten 
% Geographische Breite 10°-70°


beta = [-2.034 ; -3.7106; -4.923; -5.54]; % Abschußwinkel bzgl. x-Achse
maxcorr =  [8000 ; 8000; 8000; 8000];
for k = 1:4
    % Abschußwinkel in °
    eta = beta(k);
    phi = (k-1)*20+10;
    alpha = 0.5*asind(maxcorr(k)*g/v0^2/cosd(eta));  
    
    alphastr(k,:)= strcat("\alpha=",string(num2str(alpha,'%4.2f°')));
    betastr(k,:) = strcat("\beta=",string(num2str(eta,'%5.3f°')));
    phistr(k,:)  = strcat("\phi=",string(num2str(phi,'%2d °')));
    v0x   = v0*cosd(alpha)*cosd(eta);   % Schußrichtung Süden
    v0y   = v0*cosd(alpha)*sind(eta);   % Richtung Osten
    v0z   = v0*sind(alpha);
    % Wurfzeit und Wurfweite ohne Reibung
    tmax = 2*v0z/g;
    t = linspace(0,tmax,maxz);
    xmax(k) = v0^2*sind(2*alpha)*cosd(eta)/g;
    % Wurfzeit ohne Reibung
    t = linspace(0,tmax,maxz);
    cphi = cosd(phi);
    sphi = sind(phi);
    A   = -v0x*sphi-v0z*cphi;
    B   = v0y-g*cphi/2/omeg;
    arg = 2*omeg*t;
    % Höhe
    zw(k,:) = v0z*t -0.5*g*t.^2 + (A*(t-sin(arg)/2/omeg) - ...
              B*(1-cos(arg))/2/omeg)*cphi + 0.5*g*t.^2*cphi*cphi;
    % Ost-West-Richtung
    yw(k,:)  = (A*(1-cos(arg)) - B*(2*omeg*t-sin(arg)))/2/omeg; 
    % Süd-Nord-Richtung
    xw(k,:)   = v0x*t+(A*(t-sin(arg)/2/omeg)-B*(1-cos(arg))/2/omeg)*sphi+...
                0.5*g*t.^2*cphi*sphi;
end

alpha0 =  0.5*asind(maxcorr(k)*g/v0^2);  
alphastr(5,:)= strcat("\alpha=",string(num2str(alpha0,'%4.2f°')));
betastr(5,:) = "\beta=0°";
phistr(5,:)  = "ohne   ";
zw(5,:) = v0z*t -0.5*g*t.^2;
xw(5,:) = v0x*t; 
yw(5,:) = 0; 
 
for iPlot = 1:5
  for k=1:length(zw(iPlot,:))
   if zw(iPlot,k) < -5
      zw(iPlot,k)= NaN; 
      xw(iPlot,k)= NaN; 
      yw(iPlot,k)= NaN; 
   end
  end  
end


figure();
subplot(2,1,1);
for k=1:4 
    plot(xw(k,:)-8000,zw(k,:),'LineWidth',1,'Color',Colors(k,:),'LineStyle',Style(k));
    hold on;
end
plot(xw(5,:)-8000,zw(5,:),'LineWidth',1,'Color',Colors(9,:),'LineStyle','-');
grid on;
grid on;
xlabel('Weite in m ','FontSize',14);
ylabel('Höhe in m ','FontSize',14);
ylim([-10,10]);
xlim([-10,10]);
legend(phistr, 'location','southwest');
legend box off;
set(gca,'Fontsize', 16);

subplot(2,1,2);
for k=1:4 
    plot(xw(k,:)-8000,yw(k,:),'LineWidth',2,'Color',Colors(k,:),'LineStyle',Style(k));
    hold on;
end
plot(xw(5,:)-8000,yw(5,:),'LineWidth',2,'Color',Colors(9,:),'LineStyle',Style(5));
grid on;
grid on;
ylabel('Ostabweichung in m ','FontSize',14);
xlabel('Weite in m ','FontSize',14);
ylim([-4,inf]);
xlim([-inf,1000]);
legend(strcat(strcat(phistr, alphastr),betastr), 'location','southwest');
legend box off;
set(gca,'Fontsize', 16);



%%
% Aufgabenteil 2 c
% Optimierung Abschuß auf Ziel in Schußweite auf 8000 m
% x nach Süden y nach Osten 
% Geographische Breite 10°

phiAb = 10;
phiAbstr = strcat("\phi=",string(num2str(phiAb,'%2d °')));
beta = [-2.026 ; -2.034; -2.038; -2.046]; %Abschußwinkel bzgl. x-Achse
 
for k = 1:4
    % Abschußwinkel in °
    eta = beta(k);
    betastr(k,:) = strcat("\beta=",string(num2str(eta,'%5.3f°')));
    alpha = 0.5*asind(8000.05*g/v0^2/cosd(eta));
    alphastr = strcat("\alpha=",string(num2str(alpha,'%4.2f°')));
    v0x   = v0*cosd(alpha)*cosd(eta);   %Schußrichtung Süden
    v0y   = v0*cosd(alpha)*sind(eta);   %Richtung Osten
    v0z   = v0*sind(alpha);
    % Wurfzeit und Wurfweite ohne Reibung
    tmax = 2*v0z/g*1.1;
    t = linspace(0,tmax,maxz);
    xmax(k) = v0^2*sind(2*alpha)*cosd(eta)/g;
    % Wurfzeit ohne Reibung
    t = linspace(0,tmax,maxz);
    cphi = cosd(phiAb);
    sphi = sind(phiAb);
    A   = -v0x*sphi-v0z*cphi;
    B   = v0y-g*cphi;
    arg = 2*omeg*t;
    % Höhe
    zw(k,:) = v0z*t -0.5*g*t.^2 + (A*(t-sin(arg)/2/omeg) - ...
              B*(1-cos(arg))/4/omeg^2 + 0.5*B*t.^2)*cphi;
    % Ost-West-Richtung
    yw(k,:)  = (A*(1-cos(arg)) - B*(t-sin(arg)/2/omeg))/2/omeg; 
    % Süd-Nord-Richtung
    xw(k,:)   = v0x*t+(A*(t-sin(arg)/2/omeg)-B*(1-cos(arg))/4/omeg^2+...
                0.5*B*t.^2)*sphi;
end

betastr(5,:) = "ohne Coriolis";
zw(5,:) = v0z*t -0.5*g*t.^2;
xw(5,:) = v0x*t; 
yw(5,:) = 0; 

figure();
subplot(2,1,1);
for k=1:4 
    plot(xw(k,:)-8000,zw(k,:),'LineWidth',1,'Color',Colors(k,:),'LineStyle',Style(k));
    hold on;
end
plot(xw(5,:)-8000,zw(5,:),'LineWidth',1,'Color',Colors(9,:),'LineStyle','-');
grid on;
grid on;
xlabel('Zielweitenabweichung in m ','FontSize',14);
ylabel('Höhe in m ','FontSize',14);
title(phiAbstr+alphastr,'FontWeight','normal','FontSize',14);
ylim([-2,2]);
xlim([-10 10]);
legend(betastr, 'location','southwest');
legend box off;
set(gca,'Fontsize', 16);

subplot(2,1,2);
for k=1:4 
    plot(xw(k,:)-8000,yw(k,:),'LineWidth',2,'Color',Colors(k,:),'LineStyle',Style(k));
    hold on;
end
plot(xw(5,:)-8000,yw(5,:),'LineWidth',2,'Color',Colors(9,:),'LineStyle',Style(5));
grid on;
grid on;
title(phiAbstr+alphastr,'FontWeight','normal','FontSize',14);
ylabel('Ostabweichung in m ','FontSize',14);
xlabel('Zielweitenabweichung in m ','FontSize',14);
ylim([-0.1,0.1]);
xlim([-2,2]);
legend(betastr, 'location','southeast','NumColumns',3);
legend box off;
set(gca,'Fontsize', 16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

