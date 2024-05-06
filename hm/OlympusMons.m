% -------------------------------------------------------------------------
% OlympusMons.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial basiert z.T. auf 
% Michael Allison und Megan McEwan
% Planetary and Space Science 48(2000) 215
% -------------------------------------------------------------------------
% Sonnenaufgaenge am Olympus Mons
% 
% Programm berechnet zunaechst die Zeitgleichung auf dem Mars und daraus
% anschliessend Sonnenaufgaenge am Olympus Mons im Jahr 2030, mit und ohne
% Beruecksichtigung seiner Hoehe.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung Epoche
dt1 = datetime('2006-07-16 11:02:00');
T1  = juliandate(dt1);
t1  = Jd2JJht(T1);

Aequi = 'J2000';  % Es pielt keine Rolle ob J2000 oder Datum, da die KOS
                  % fuer wahre und mittlere gleich sind und nur die
                  % Differenz eine Rolle spielt.
[BaPa,BaPadot]= OrbitParameter(T1,Aequi); 
EarthDay      = 86400; % Laenge eines Erdtages 
 
%%
%%-----------------------------------------------------------------------%%
% Beginn Rechnung

PlNr=4; % Mars, es kann auch ein anderer Planet genommen werden

% Bahnparameter aus Datenfiles
ex    = BaPa.eP(PlNr);          % Exzentrizitaet
exdot = BaPadot.edot(PlNr);     % Exzentrizitaetsveränderung
epsP  = BaPa.epsP(PlNr);        % Achsneigung
M0P   = BaPa.M0P(PlNr);         % mittl. Anaomalie bei J2000
n_anom= BaPa.NP(PlNr);          % anomal. Umlauf
n_trop= BaPadot.L0dot(PlNr);    % trop. Umlauf
M0dot = BaPadot.M0dot(PlNr);    % trop. Umlauf 
L0P   = BaPa.L0P(PlNr);         % mittl. Laenge bei J2000
PiqP  = BaPa.PiqP(PlNr);        % Laenge des Perihels bei J2000 
Piqdot= BaPadot.Piqdot(PlNr);   % Praezession des Perihels 

% Werte müssen aus der Literatur beschafft werden oder berechnet werden
LSP      = 250.9999;            % Planetozentr. Laenge (!!!) des Perihel
                                % bezogen auf Fruehlingspunkt des Planeten%
                                % Für die Erde ist einfach LSP = \omega +180
% Daraus erhaelt man:                                
alphaMS0 = LSP + M0P;          
% Faktoren fuer Mittelpunktsgleichung aus Exzentrizität
% ex     = ex + exdot*t1;
fac(1) = (2*ex^1 - 1/4*ex^3   + 5/96*ex^5);
fac(2) = (5/4*ex^2 - 11/24*ex^4 + 17/192*ex^6);
fac(3) = (13/12*ex^3 - 43/63*ex^5);
fac(4) = (103/96*ex^4 - 451/480*ex^6);
fac(5) = (1097/960*ex^5);

% Länge eines Orbits/Jahres in Julianischen Tagen 
NrJDpOrbit = 365.25*36000/n_trop;           % Erdtage eines Jahres 

% Länge eines Sonnentages (sol) und Länge eines Orbits/Jahres in sols
SidTag        = BaPa.SidTagP(PlNr);         % Siderischer Tag
SolTag        = SidTag/(1-SidTag/(NrJDpOrbit*86400));    % Solartag auf Planeten (sol)
NrSolDpOrbit  = NrJDpOrbit*86400/SolTag;    % Solartage eines Jahres 
SolSidRatio   = SolTag/SidTag;              % Verhältnis Sol/Sid Tag Mars
SolJDRatio    = SolTag/EarthDay;            % Verhältnis Sol/Erdtag
SidJDRatio    = SidTag/EarthDay;            % Verhältnis Sid/Erdtag

% Datumsvektor für Zeitbereich
NrDays      = floor(NrSolDpOrbit)+1;              % Sols für Berechnung
Tv          = linspace(T1,T1+NrJDpOrbit+1,NrDays+1);  % Jul. Tage
tv          = Jd2JJht(Tv);                      % Julian. Jahrhunderte 
dtv         = datetime(Tv,'convertfrom','juliandate');
index = 1;

for tag = 1:length(dtv)  % in Julianischen Tagen
    MonatsTag = day(dtv(tag),'dayofmonth');
    if MonatsTag == 1
      BeginnMonat(index) = tag;
      Monat(index) = month(dtv(tag));
      Jahr(index)  = year(dtv(tag))-2000;
      index = index+1;
    end
end

% Schleife ueber die Julian. Tage eines Jahres
for tag = 1:NrDays+1  % in Julianischen Tagen
    tz = tv(tag);
    M(tag)       = M0P + n_anom*tz;             % Mittl. Anomalie
    alphaMS(tag) = alphaMS0 + n_trop*tz;        % RA der mittl. Sonne
    MPG(tag) = 0;
    for nz=1:5                                  % MPG
        MPG(tag)  = MPG(tag)+ rad2deg(fac(nz)*sind(nz*M(tag)));
    end
    LS(tag)  = alphaMS(tag) + MPG(tag);         % Planetozentr. Laenge
    R2E(tag) = 0;
    for nz=1:3                                  % Reduktian auf Aequator
        faktor   = (-(tand(epsP/2))^2)^nz;
        R2E(tag) = R2E(tag)+rad2deg(faktor*sind(2*nz*LS(tag))/nz);
    end
    ZGL(tag) = -4*(MPG(tag)+ R2E(tag));         % ZGL
    Dec(tag) = asind(sind(epsP)*sind(LS(tag)))+ 0.25*sind(LS(tag));  % Dekl.
end



%%
% Bild Zeitgleichung
figure(); 
plot(dtv,ZGL,...
            'Color', Colors(4,:),'LineStyle','-','LineWidth',2);
ylim([-60 +60]);
ylabel('Datum  UT')
ylabel('ZGL in min');
grid on;
set(gca,'Fontsize',16);

% Bild Analemma
figure();
% plot(ZGL,Dec,...
%         'Color', Colors(4,:),'LineStyle','-','LineWidth',2);
plot(ZGL, Dec,'-+','MarkerIndices',BeginnMonat,'LineWidth',2,...
      'Color',Colors(4,:),'LineStyle','-');
for k=1:length(BeginnMonat)
        mylabels(k,:)=sprintf('1.%02u.%02u',Monat(k),Jahr(k));
        xL(k)=ZGL(BeginnMonat(k));
        yL(k)=Dec(BeginnMonat(k));
        LabelPoints(xL(k), yL(k) ,mylabels(k,:),'e',0.3,0,'FontSize',...
                    12,'Color',Colors(4,:));
end

grid on,
ylim([-30 30]);
xlim([-60 60]);
ylabel('Deklination')
xlabel('ZGL in min');
set(gca,'Fontsize',16);


%%
% Initialisierung fuer Sonnenaufgang

% Koordinaten Olympus Mons
lambdaOM=-133.8;
phiOM=18.65;
heightOM = 21.2; % Hoehe in km
rM = 3389.5; % Mars-Radius in km 
lambda = lambdaOM;
phi = phiOM;

% Marszeit
sid = SidJDRatio; % Siderischer Tag in Tagen
sol = SolJDRatio; % sol in Julian. Tagen

sinphi=sind(phi);
cosphi=cosd(phi);

% Sonnendurchmesser vom Mars 19.3'
% Parallaxe 3"
% vernachlaessige Refraktion in Marsatmosphaere
hr = -9.6/60; %hr Korrektur 

% Arraygroessen festlegen
zeitdiff=zeros(1,670);
taglaenge=zeros(1,670);
auf=zeros(1,670);
auf_top=zeros(1,670);
aufJD=zeros(1,670);
unter=zeros(1,670);
q =zeros(1,670);


% Beginn Rechnung Sonnenaufgaenge-------------------------------------------

% Zeitvektor mit Abstand 1 sol und Laenge 1 Mars-Jahr in JD
T_v = linspace(T1,T1+669*sol,670);
% Berechnung Kulminationen in MOZ aus Zeitgleichung
kulmi = 12 + ZGL/60;


% Berechnung Sonnenaufgang
q=(sind(hr)-sind(phi)*sind(Dec))./(cosd(phi)*cosd(Dec));

% Gipfel des Olympus Mons
el_corr = acosd(rM/(rM+heightOM)); % Hoehenkorrektur
disp(el_corr/15);
kulmi_top = kulmi - el_corr/15;

for k = 1:length(T_v)
    zeitdiff(k)    = acosd(q(k))/15; 
    auf(k)    = kulmi(k) - zeitdiff(k);
    auf_top(k) = kulmi_top(k) -zeitdiff(k);
    unter(k)  = kulmi(k) + zeitdiff(k);
    taglaenge(k)   = 2*zeitdiff(k);
end

% Graphische Ausgabe--------------------------------------------------------

figure();
XTime = (T_v-T_v(1))/sol;

plot(XTime,auf,'Color',Colors(2,:),'LineWidth',2);
hold on;
plot(XTime,auf_top,'Color',Colors(2,:),'LineWidth',2,'LineStyle','--');
hold on;

txt = "Sonnenaufgang am Fuss und Gipfel des Olympus Mons";
text(20,25, txt,'FontSize',18);

grid on;
xlim([0,668]);
lgd=xlabel('Zeit in sols (Marstagen)');
lgd.FontSize=16;

ylim([4,8]);
lgd=ylabel('Zeit in Lokalzeit (MOZ)');
lgd.FontSize=16;
yticks manual;
yticks([4 5 6 7 8]);
yticklabels({'04:00' '05:00' '06:00' '07:00' '08:00'});
lgd=legend('SA (Fuss)','SA (Gipfel)',...
    'location','southeast');
legend boxoff;
lgd.FontSize=16;
hold on;
set(gca,'Fontsize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
