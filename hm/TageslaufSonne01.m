% -------------------------------------------------------------------------
% TageslaufSonne01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Tageslauf der Sonne 
% Der Tageslauf der Sonne wird nach der Näherungsformel (3.84) berechnet. 
% Die geographische Breite kann variiert werden. Es wird die Höhe der Sonne 
% für verschiedene Jahreszeiten und geografische Breiten dargestellt. 
% Die ZGL wird nicht berückcksichtigt.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
t1 = datetime('2000-01-01 0:0:00');
% Julianisches und Modifiziertes Julianisches Datum
T = juliandate(t1);
MJuDa =juliandate(t1,'modifiedjuliandate');

%Schiefe der Ekliptik J2000;
eps=deg2rad(23.43929111); 
%Refraktion und scheinbarer Durchmesser in °
hr=50/60;

%Äquinoktien und Sonnenwenden
tag(1)=datetime(2000,3,21);
tag(2)=datetime(2000,6,20);
tag(3)=datetime(2000,9,22);
tag(4)=datetime(2000,12,21);

%Vorbereitung Felder
hoehe1    = zeros(4,240);
hoehe2    = zeros(4,240);
hoehe3    = zeros(4,240);
hoehe4    = zeros(4,240);
horizont  = zeros(240);
titlestr1 = strings(4,1);
titlestr2 = strings(4,1);
tz        = zeros(1,4);
zeit      = zeros(1,240);

%_______________________________________________________
%Beginn Hauptprogramm


% Schleife über vier Breitengrade
for k=1:4
  phi=(k-1)*25;
  % Schleife über vier Zeitpunkte
  for m=1:4
    tz(m)=day(tag(m),'dayofyear');
    % Berechnung Deklination mit Sinus-Näherungsformel für Deklination
    % festgehalten täglich bei 12:00 UT
    for n=1:240
        tau =((n-1)-120)*2*pi/240;
        zeit(n) = (n-1)/10;
        delta = rad2deg(0.4095)*sin(0.016906*((tz(m)-12/24+zeit(n)/24)-80.086));
        switch k
        case 1 
           hoehe1(m,n)=asind(sind(phi)*sind(delta)+cos(tau)*cosd(phi)*cosd(delta));
        case 2 
           hoehe2(m,n)=asind(sind(phi)*sind(delta)+cos(tau)*cosd(phi)*cosd(delta));
        case 3
           hoehe3(m,n)=asind(sind(phi)*sind(delta)+cos(tau)*cosd(phi)*cosd(delta));
        case 4 
           hoehe4(m,n)=asind(sind(phi)*sind(delta)+cos(tau)*cosd(phi)*cosd(delta));
        end   
    end
    titlestr1(m)=datestr(tag(m),' dd.mm.');
  end
end

%Ausgabe
figure('Name','Taglauf der Sonne für verschiedene Jahreszeiten');

for k=1:4
    subplot(2,2,k);
    hold on;
    switch k
        case 1 
            for m=1:4
                plot(zeit,hoehe1(m,:),'Color', Colors(m+1,:),'Linewidth',2);
             end
         case 2
             for m=1:4
                plot(zeit,hoehe2(m,:),'Color', Colors(m+1,:),'Linewidth',2);
             end
         case 3 
             for m=1:4
                plot(zeit,hoehe3(m,:),'Color', Colors(m+1,:),'Linewidth',2);
             end
         case 4
             for m=1:4
                plot(zeit,hoehe4(m,:),'Color', Colors(m+1,:),'Linewidth',2);
             end
    end
    plot(zeit,horizont,'k--'); 
    grid on;
    phi=(k-1)*25;
    phistr =sprintf(' = %+4.1f°',phi);
    ttl=title(['Breite {\phi}',phistr]); 
    lgd=legend(titlestr1(1),titlestr1(2),titlestr1(3),titlestr1(4),'location','south');
    legend boxoff;
    xlim([1 24]);
    ylim([-90 90]);
    xticks(0:2:24);
    yticks(-90:30:90);
    xlabel('Zeit in h');
    ylabel('Hoehe in °');
    lgd.FontSize = 12;
    ttl.FontSize = 12;
    set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
