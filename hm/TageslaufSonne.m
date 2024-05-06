% -------------------------------------------------------------------------
% TageslaufSonne.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Tageslauf der Sonne 
% Berechnung des Tagbogenverlaufs der Sonne über dem Horizont.
% Es wird einmal die Deklination der Sonne nach einer Näherungsformel 
% berechnet und zum anderen "numerisch" unter Berücksichtigung der 
% Refraktion und Zeitgleichung mit der Kepler-MPG-Lösung.
% Die geographische Breite kann variiert werden.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
t1 = datetime('2000-01-1 0:0:00');
% Julianisches und Modififziertes Julianisches Datum
T = juliandate(t1);
MJuDa =juliandate(t1,'modifiedjuliandate');

%Schiefe der Ekliptik J2000;
eps=deg2rad(23.43929111); 
%Refraktion und scheinbbarer Durchmesser in 0
hr=50/60;
%Geograph. Breite in 0
phi=50;
phistr =sprintf(' = %+4.1f0',phi);

%Sonnenwenden und Tagundnachtgleichen
tag(1)=datetime(2000,3,20);
tag(2)=datetime(2000,6,20);
tag(3)=datetime(2000,9,22);
tag(4)=datetime(2000,12,21);

%Vorbereitung Felder
hoehe1    = zeros(4,240);
hoehe2    = zeros(4,240);
titlestr1 = strings(4,25);
titlestr2 = strings(4,25);

%Beginn Hauptprogramm_______________________________________________________
% 
% Schleife ueber 4 Zeitpunkte des Jahres
MJuDa0=MJuDa-1;
T0=T-1;

zeit = zeros(1,240);
for m=1:4
    tz(m)=day(tag(m),'dayofyear');
    % Berechnung Deklination mit Sinus-Näherungsformel für Deklination
    % festgehalten für einen Tag bei 12:00 UT
    delta = rad2deg(0.4095)*sin(0.016906*(tz(m)-80.086));
    for n=1:240
        tau =((n-1)-120)*2*pi/240;
        zeit(n) = (n-1)/10;
        delta = rad2deg(0.4095)*sin(0.016906*((tz(m)-12/24+zeit(n)/24)-80.086));
        %KOS Trafo
        hoehe1(m,n)=asind(sind(phi)*sind(delta)+cos(tau)*cosd(phi)*cosd(delta));
    end
    titlestr1(m)=sprintf(' Näherung');
end

% Schleife ueber 4 Zeitpunkte
for m=1:4
  MJuDa=MJuDa0+day(tag(m),'dayofyear');
  T=T0+day(tag(m),'dayofyear');
% Berechnung Deklination numerisch nach MPG aller 6 min
  for n=1:240
      zeit(n) = (n-1)/10;
      MJD=MJuDa+zeit(n)/24;
      TJD=T+zeit(n)/24;
      %Berechnung numerisch nach MPG 
      [RA,Dec]=KeplerSonne(TJD,eps); 
      %Berechnung des Stundenwnkels
      tau=GMST(MJD)-RA;
      tau = wrapToPi(tau);
      %KOS Trafo
      sinusalt = sind(phi)*sin(Dec)+cosd(phi)*cos(Dec)*cos(tau);   
      h=asind(sinusalt);
      %Refraktionskorrektur
      hoehe2(m,n)=h+1.02./tand(h+10.3./(h+5.11))/60;
      if hoehe2(m,n) < -1 
          hoehe2(m,n)= NaN;
      end
  end
  titlestr2(m)=strcat(datestr(tag(m),' dd.mm.2000 '), ' exakt');
end



%Ausgabe
figure('Name','Taglauf der Sonne Naeherung und exakt)');
subplot(1,2,1);
hold on;

for k=1:4
    plot(zeit,hoehe2(k,:),'Color',Colors(k+1,:),'LineWidth',1);
end
for k=1:4
    plot(zeit,hoehe1(k,:),':','Color',Colors(k+1,:),'LineWidth',1);
end
ttl=title(['Taglauf der Sonne bei {\phi}',phistr]); 
ttl.FontSize = 14;
xlim([1 24]);
xticks(0:2:24);
lgd=xlabel('Zeit in h');
lgd.FontSize = 14;
lgd=ylabel('Hoehe in °');
lgd.FontSize = 14;
lgd=legend(titlestr2(1),titlestr2(2),titlestr2(3),titlestr2(4),titlestr1(1),titlestr1(2),titlestr1(3),titlestr1(4),'Location','northwest','NumColumns',2);
lgd.FontSize = 12;
legend boxoff;
grid on;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

subplot(1,2,2);
hold on;
for k=1:4
    plot(zeit,hoehe2(k,:)-hoehe1(k,:),':','Color',Colors(k+1,:),'LineWidth',1);
    titlestr2(k)= datestr(tag(k),' dd.mm.2000 ');
 end
ttl=title('Abweichung der Näherung in °');
ttl.FontSize = 14;
xlim([1 24]);
xticks(0:2:24);
lgd=xlabel('Zeit in h');
lgd.FontSize = 14;
lgd=ylabel('Abweichung Hoehe in °');
lgd.FontSize = 14;
lgd=legend(titlestr2(1),titlestr2(2),titlestr2(3),titlestr2(4));
lgd.FontSize = 12;
legend boxoff;
grid on;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------