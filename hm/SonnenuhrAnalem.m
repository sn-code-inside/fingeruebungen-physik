% -------------------------------------------------------------------------
% Wandsonnenuhr.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Lineatur einer analemmatischen Sonnenuhr
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Position von Görlitz Germany
lati  = 51.2;    % Geograph. Breite
longi = 15;    % Geograph. Länge
latlong = sprintf(' Breite %4.1f° Länge %4.1f° (Uhrzeit in UT)', lati, longi);
phi = deg2rad(lati);
maanden = [31,28,31,30,31,30,31,31,30,31,30,31]; % months

%Schiefe der Ekliptik J2000;
eps=deg2rad(23.43929111); 


%Vorbereitung Felder
hoehe    = zeros(6,24);
horizont  = zeros(1,48);

hend = 13;
dend = 6;
titlestr1 = strings(dend,1);
titlestr2 = strings(hend,1);
titlestr3 = strings(hend+dend,1);
tz        = zeros(1,6);
zeit      = zeros(1,24);
%%ausgewählte Tage des Jahres
tag(1)=datetime(2000,1,1);
tag(2)=datetime(2020,3,1);
tag(3)=datetime(2020,5,1);
tag(4)=datetime(2020,7,1);
tag(5)=datetime(2020,9,1);
tag(6)=datetime(2020,11,1);
H = 3;
%_______________________________________________________
%Beginn Hauptprogramm

% Schleife über 6 Tage und 24 Zeitpunkte
for h=1:48 %Stunden des Tages
    MOZ(h) = (h+2)/2;  %MOZ für Görlitz
    for n=1:6  %ausgewählte Tage des Jahres
        %Berechnet Alt/Azimuth
        [az, alt] = azialti(2020, (n-1)*2+1, 1, MOZ(h)-1, 0, lati, longi); 
        az0h(h, n) = az;
        alt0h(h, n)= alt;
        titlestr1(n)=datestr(tag(n),' dd.mm.');
    end
end

for h=1:48 %Stunden des Tages
    for n=1:6  %ausgewählte Tage des Jahres
        %Berechnet Alt/Azimuth
        fac = H*cotd(alt0h(h, n));
        xns1 (h,n) = fac*cosd(az0h(h,n));
        yns1 (h,n) = fac*sind(az0h(h,n));
    end
end
for h=1:48
    for n=1:6
         if alt0h(h,n) <= 0 
           yns1(h,n) = NaN;  
           xns1(h,n) = NaN;  
         end
    end
end

%Ausgabe Linien konstanter Deklination
figure('Name',' konstanter Deklination');
hold on;
for h=1:48
for n=1:dend
        pl1(n)= plot(yns1(h,n),xns1(h,n),'+','Color',...
                Colors(n+1,:),'Linewidth',2);
end
end
for h=10:2:36
    text(yns1(h,4)-0.2,xns1(h,4)-0.5,strcat(num2str(MOZ(h),2),'h'),'Color',...
                Colors(4+1,:));
    text(yns1(h,1)-0.2,xns1(h,1)+0.5,strcat(num2str(MOZ(h),2),'h'),'Color',...
                Colors(1+1,:));
end
for n=1:dend
        plot(yns1(:,n),xns1(:,n),'Color', Colors(n+1,:),...
            'Linewidth',2,'LineStyle',':');
end
grid on;
PlotCircle (0,0,0.25,'k',2);  %Nodus
phistr =sprintf(' = %+4.1f°',lati);
ttl=title('Görlitz:  Länge 15° O, Breite 51.2° N'); 
lgd=legend(pl1,titlestr1,'location','bestoutside','NumColumns',1);
legend boxoff;
axis equal
ylim([-5 15]);
xlim([-10 10]);
xlabel('WO in m');
ylabel('NS in m');
lgd.FontSize = 12;
ttl.FontSize = 12;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% Schleife über 365 Tage und 7 Zeitpunkte
for h=1:hend %Stunden des Tages
    for n=1:365  %Tage des Jahres
        d = n;
        m = 1;
        while (d>maanden(m))  % Fortlaufende Tage des Jahres 
                              % umgerechnet in Tage des Monats
            d = d - maanden(m); % Monatstag
            m = m+1;            % Monat
        end
        GMT(h)=MOZ(h*2)+4-1;
        %Berechnet Alt/Azimuth
        [az, alt] = azialti(2020, m, d, GMT(h), 0, lati, longi); 
        azA(h, n) = az;
        altA(h, n)= alt;
    end
    titlestr2(h)= strcat('MOZ : ',num2str(GMT(h)+1));
    titlestr2(h)= strcat(titlestr2(h), ':00 h');
end

for h=1:hend %Stunden des Tages
    for n=1:365 %ausgewählte Tage des Jahres
        %Berechnet Alt/Azimuth
        fac = H*cotd(altA(h, n));
        xns2 (h,n) = fac*cosd(azA(h,n));
        yns2 (h,n) = fac*sind(azA(h,n));
    end
end
for h=1:hend
    for n=1:365
         if altA(h,n) < 0 
           yns2(h,n) = NaN;  
           xns2(h,n) = NaN;  
         end
    end
end

%Ausgabe der Stundenlinien (Analemma)
figure1=figure('Name','Stundenlinien (Analemma)');
hold on;
SSW =172; %Sommersonnenwende
WSW =355; %Wintersonnenwende

for h=1:hend
   pl1(h)= plot(yns2(h,:),xns2(h,:),...
             'color',Colors(h,:),'Linewidth',1);
   plot(yns2(h,SSW),xns2(h,SSW),'o',...
             'color',Colors(h,:),'Linewidth',1);
   plot(yns2(h,WSW),xns2(h,WSW),'d',...
             'color',Colors(h,:),'Linewidth',1);
   text(yns2(h,SSW)-0.2,xns2(h,SSW)-0.5,strcat(num2str(GMT(h)+1,2),'h'),...
             'Color', Colors(h,:));
end
grid on;
phistr =sprintf(' = %+4.1f°',lati);
ttl=title('Görlitz:  Länge 15° O, Breite 51.2° N'); 
PlotCircle (0,0,0.25,'k',2);  %Nodus
lgd=legend(pl1,titlestr2,'location','bestoutside','NumColumns',1);
legend boxoff;
axis equal
ylim([-5 15]);
xlim([-10 10]);
xlabel('WO in m');
ylabel('NS in m');
lgd.FontSize = 12;
ttl.FontSize = 12;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

%Ausgabe Lineatur
figure1=figure('Name','Lineatur');
axis equal
ylim([-5 25]);
xlim([-15 15]);
hold on
for n=1:dend
    pl3(n)= plot(yns1(:,n),xns1(:,n),'+','Color', Colors(n+1,:),...
            'Linewidth',1);
    titlestr3(n)=titlestr1(n);
    plot(yns1(:,n),xns1(:,n),'Color', Colors(n+1,:),'Linewidth',1,...
        'Linestyle',':');
end
for h=10:2:36
    text(yns1(h,4)-0.2,xns1(h,4)-0.5,strcat(num2str(MOZ(h),2),'h'),'Color',...
                Colors(4+1,:));
    text(yns1(h,1)-0.2,xns1(h,1)+0.5,strcat(num2str(MOZ(h),2),'h'),'Color',...
                Colors(1+1,:));
end
for h=1:hend
    pl3(n+h)= plot(yns2(h,:),xns2(h,:),...
              'color',Colors(h,:),'Linewidth',1);
    plot(yns2(h,SSW),xns2(h,SSW),'o',...
             'color',Colors(h,:),'Linewidth',1);
    plot(yns2(h,WSW),xns2(h,WSW),'d',...
             'color',Colors(h,:),'Linewidth',1);
    titlestr3(n+h)=titlestr2(h);
end   
grid on;
phistr =sprintf(' = %+4.1f°',lati);
ttl=title('Görlitz:  Länge 15° O, Breite 51.2° N'); 
PlotCircle (0,0,0.25,'k',2);  %Nodus
lgd=legend(pl3,titlestr3,'location','bestoutside','NumColumns',1);
legend boxoff;
axis equal
ylim([-5 20]);
xlim([-15 15]);
xlabel('WO in m');
ylabel('NS in m');
lgd.FontSize = 12;
ttl.FontSize = 12;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Funktionen
% -------------------------------------------------------------------------


% Berechnung von Höhe und Azimuth
function [azi,alti] = azialti(y, m, d, h, mins, lati, longi) 
    h = h+mins/60.0;

    % Part I: Inspired by http://www.stargazing.net/kepler/sun.html
    %1. Find the days before J2000.0 (d) 
       d =  367 * y - floor(7 * (y +floor( (m + 9)/12))/ 4) + floor(275 * m / 9) + d - 730531.5 + h / 24;
   
    %2. Find the Mean Longitude (L) of the Sun (siehe Formel 3.65 im Buch) 
       L = wrapTo360(280.461 + 0.9856474 * d);
      
    %3. Find the Mean anomaly (g) of the Sun
       g = 357.528 + 0.9856003 * d;
      
    %4. Find the ecliptic longitude (lambda) of the sun
        % (Mittelpunktsgleichung)
       lambda = L + 1.915 * sind(g) + 0.020 * sind(2*g);
   
    %5. Find the obliquity of the ecliptic plane (epsilon)
       epsilon = 23.439 - 0.0000004 * d;
      
    %6. Find the Right Ascension (alpha) and Declination (delta) of Sun
        % Formel 3.4 im Buch
       Y = cosd(epsilon) * sind(lambda);
       X = cosd(lambda);
       %Rektaszension    
       alpha = rad2deg(atan2(Y, X));     
       %Deklintaion
       delta = asind(sind(epsilon)*sind(lambda));
      
   % Part II: From http://www.geoastro.de/elevaz/basics/index.htm
   %compute Sidereal time at Greenwich:
   % (Naherungsformel)
       T = d/36525.0;  % Zeit in Julianischen Jahrhunderten
       theta0 = wrapTo360(280.46061837 + 360.98564736629*d + 0.000387933*T*T - T*T*T/38710000.0);
       theta = theta0 + longi;
       tau = theta - alpha; %Stundenwinkel)
       beta = lati;
       % Formel 3.4 im Buch
       alt = asin(sind(beta)*sind(delta) + cosd(beta)*cosd(delta)*cosd(tau));
       az  =  atan2(- sind(tau) , (cosd(beta)*tand(delta) - sind(beta)*cosd(tau)));
       alti= rad2deg(alt);
       azi = wrapTo180(rad2deg(az)-180);
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
