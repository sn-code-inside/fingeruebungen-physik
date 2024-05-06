% -------------------------------------------------------------------------
% TageslaufMond.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Tageslauf des Mondes 
% Berechnung des Tagbogenverlaufs des Mondes ueber dem Horizont.
% Es wird die Deklination des Mondes 
% "numerisch" unter Beruecksichtigung der Parallaxe und
% Refraktion mit der Kepler-MPG-Lösung berechnet.
% Die geographische Breite kann variiert werden.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
dt1    = datetime('today');
dt1    = datetime('2020-03-01 00:00:00');
dt1    = datetime('2018-03-01 00:00:00');
% Julianisches und Modififziertes Julianisches Datum
T1    = juliandate(dt1);
MJuDa = juliandate(dt1,'modifiedjuliandate');
MJuDa0=MJuDa-1;
T0=T1-1;

%Schiefe der Ekliptik J2000;
epsE    = EpsErde(T1);  
epsErad = deg2rad(epsE); 

%Refraktion,Parallaxe und scheinbarer Durchmesser in °
hr = -8/60;
%Geograph. Breite/Laenge in °
phi    = 47.5;
lambda = 20;
phistr = sprintf(' = %+4.1f°',phi);
lamstr = sprintf(' = %+4.1f°',lambda);
Zone   = 2;


%%
%Beginn Hauptprogramm
% 
% Schleife ueber Zeitpunkte 
Tage =26;
NrJDpOrbit = 27.32166;
for k=1:Tage
  MJuDa=MJuDa0+k;
  T1(k)=T0+k;
  plotstr(k)=string(datetime(T1(k),'ConvertFrom','juliandate'),'dd.MM.yy');
% Berechnung Deklination numerisch nach MPG jede min
  for n=1:24*60*6
      uhrzeit(n) = (n-1)/24/60/6;
      MJD        = MJuDa + uhrzeit(n);
      TJD(k,n)   = T1(k) + uhrzeit(n);
      %Berechnung numerisch in verschiedenen Näherunge 
      Moon    = KeplerMond(TJD(k,n),epsErad); %Mondkoordinaten in Näherung
      MoonM   = KeplerMondM(TJD(k,n),epsErad); %Mondkoordinaten in Näherung
      RA     = Moon.equ(2);             %Geozentr-aequatorial
      Dec    = Moon.equ(3);             %Geozentr-aequatorial
      RAM    = MoonM.equ(2);            %Geozentr-aequatorial
      DecM   = MoonM.equ(3);            %Geozentr-aequatorial
      %Berechnung des Stundenwinkels       
      tau  = GMST(MJD)-RA+deg2rad(lambda);  
      tau  = wrapTo360(rad2deg(tau));
      tauM = GMST(MJD)-RAM+deg2rad(lambda);  
      tauM = wrapTo360(rad2deg(tauM));
      %KOS Trafo
      sinusalt = sind(phi)*sin(Dec)+cosd(phi)*cos(Dec)*cosd(tau);   
      h        = asind(sinusalt);
      Az(k,n)  = 180+atan2d((cos(Dec)*sind(tau)), ...
                (cos(Dec)*cosd(tau)*sind(phi)-sin(Dec)*cosd(phi)));
      Alt(k,n) = h;
      DecW(k,n) = Dec;
      sinusalt = sind(phi)*sin(DecM)+cosd(phi)*cos(DecM)*cosd(tauM);   
      hM       = asind(sinusalt);
      AzM(k,n) = 180+atan2d((cos(DecM)*sind(tauM)), ...
                (cos(DecM)*cosd(tauM)*sind(phi)-sin(DecM)*cosd(phi)));
      AltM(k,n)= hM;
      ZGL2(k,n) = (tau - tauM)*4;
%       %Probeausdruck
%       zeitstr = string(datetime(TJD,'ConvertFrom','juliandate'),...
%                'dd.MM.yy HH:mm:ss');    
%       fprintf('%s  %5.2f°  %6.2f°  %6.2f°    %7.2f°  %6.2f° \n', ...
%           zeitstr, rad2deg(RA), rad2deg(Dec), tau, Az(k,n), Alt(k,n) ); 
  end
end


%Bestimmung Kulminationszeitpunkte
for k=1:Tage
    [Altkum(k),Tkum(k)]   = max(Alt(k,:));
    [AltkumM(k),TkumM(k)] = max(AltM(k,:));
    Deckum(k) = DecW(k,Tkum(k));
    dtkum(k)  = datetime(TJD(k,Tkum(k)),'ConvertFrom','juliandate');
    dtkumM(k) = datetime(TJD(k,TkumM(k)),'ConvertFrom','juliandate');
    D(k)  = duration(string(dtkum(k),'HH:mm:ss'));
    DM(k) = duration(string(dtkumM(k),'HH:mm:ss'));
    Abst(1) = duration('00:00:00');
    if k > 1 
        Abst(k) = DM(k)-DM(k-1);
    end
    if Tkum(k) == 1
        ZGL(k) = duration('00:59:00');
    else
        ZGL(k) = -(D(k)-DM(k));
    end
    tagstr = string(datetime(TJD(k,1),'ConvertFrom','juliandate'),...
                'dd.MM.yyyy');
    fprintf('Kulminationszeit: %s %s UTC Hoehe: %5.2f° Mittl. Mond %s UTC Hoehe: %5.2f°  ZGL: %s   Abstand: %s \n',...
                tagstr, string(D(k)),Altkum(k),string(DM(k)),AltkumM(k), string(ZGL(k)), string(Abst(k)) ); 
end

figure('Name','ZGL');
hold on;
ttl=title(['ZGL Mond bei {\phi}',phistr,' {\lambda}', lamstr]); 
ttl.FontSize = 14;
plot(ZGL,Altkum,'+','LineWidth',2,'Color',Colors(6,:));
% plot(ZGL,rad2deg(Deckum),'+','LineWidth',2,'Color',Colors(6,:));
xtickformat('hh:mm');
xlim([duration('-00:45:00') duration('00:45:00')]);
xlabel('Zeit in min');
ylabel('Hoehe in °');
legend(string(dt1,'dd-MMM-yy'),'location','southeast');
legend box off;
grid on;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

figure('Name','Kulminationszeitpunkte ');
hold on;
ttl=title(['Kulminationszeitpunkte Mond bei {\phi}',phistr]); 
ttl.FontSize = 14;
dt = datetime(T1,'ConvertFrom','juliandate');
plot(dt,D,'o','LineWidth',2,'Color',Colors(6,:));
plot(dt,DM,'+','Color',Colors(3,:));
xtickformat('dd.MM.yy');
ytickformat('hh:mm');
ylim([duration('0:0:0') duration('24:0:0')]);
ylabel('Zeit in h');
xlabel('Tage');
grid on;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);


%Ausgabe
figure('Name','Tageslauf Mond Alt-Uhrzeit');
hold on;
for k=1:Tage
  p = plot(uhrzeit(:)*24,Alt(k,:),'Color',Colors(mod(k,10)+1,:),'LineWidth',1);
  if k < 11 
       p.LineStyle = '-';
  else if k < 21
       p.LineStyle = '-.';
      else
       p.LineStyle = ':';
      end
  end
end
ttl=title(['Tageslauf Mond bei {\phi}',phistr ,' {\lambda}', lamstr]); 
ttl.FontSize = 14;
xlim([0 24]);
ylim([-80 80]);
xticks(0:2:24);
xlabel('Zeit in h');
ylabel('Hoehe in °');
lgd=legend(plotstr,...
    'Location','northwest','NumColumns',10);
lgd.FontSize = 12;
legend boxoff;
grid on;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 14, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
