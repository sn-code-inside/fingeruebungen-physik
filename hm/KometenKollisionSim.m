%%-------------------------------------------------------------------------
% SonnensystemSimulationLangeZeiten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Himmelsmechanik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Simulation der Bewegung der Planeten und der Sonne 
% um den Schwerpunkts des Sonnensystem (Solar System Barycenter - SSBC)
% 
% Vergleich der numerischen Ergebnisse mit den genauen Reihen-Lösungen
% aus der Störungstheorie.
% 
% -------------------------------------------------------------------------

%% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Style = ["-", "-.", ":", "--", ":"];
Colors=GetColorLines;

% Anzahl der Himmelskörper

PNamen = ["Merkur" "Venus" "Erde" "Mars" "Jupiter" "Saturn" "Uranus" ...
          "Neptun" ... %"9P/Tempel 1" 
          "Swift-Tuttle" "Sonne"];
NK = length(PNamen);
PlanetsOnly = 8;  % Anzahl der Planeten

% Parameter, Rechnungen erfolgen alle in m/s
AE   = 1.495978707e11;            % AE in m
AEk  = 1.495978707e08;            % AE in km
G    = 6.6743e-11;                % Gravitationskonstante in m^3 / kg / s^2
Tag  = 24 * 60 * 60;              % in s
RSun = 696342e3; rsun = RSun/AE;  % Sonnenradius in m bzw AE
RE   = 6378000;  rE   = RSun/AE;  % Erdradius in m bzw AE
GAE  = 1.48600650331318e-34;      % G in AE^3/d^2/kg
MS   = 1.98840987129736e+30;      % in kg
muG  = GAE*MS;                    % µG in AE^3/d^2/kg

%% Einladen der JPL Daten Erde und Mars für 2022

% Daten in den csv-Files sind in km und km/s
% Umrechnung in m und m/s
JPLD = JPLVectorTable('JPL_Earth_SSB_Long.csv');
NSteps = length(JPLD.T);  % In Schritten von 1h
NSteps = 40000;
for k= 1:NSteps
    JPLE.r(k) = sqrt((JPLD.x(k))^2+(JPLD.y(k))^2+ (JPLD.z(k))^2)*1000;
    JPLE.xyz(:,k)=[JPLD.x(k); JPLD.y(k); JPLD.z(k)]*1000;
end

%% Vorbereitung numerische Berechnung
% Simulationszeit und Zeitschrittweite [in Tagen].
dt1str= "01.01.2022";
dt1   = datetime(dt1str);    
T1 = juliandate(dt1); % Julianisches Datum
dt1str= datestr(datetime(T1,'convertfrom','juliandate'),'mmm-yyyy');

dt    = 86400;             %Schrittweite in s
tend  = NSteps*dt;         %Ende der Simulation in s
T2  = T1+(tend-dt)/Tag;  %Ende der Simulation in Julian. Datum
t_vector=0:dt:tend-dt;
T_vector=T1:dt/Tag:T2;
dt2str = datestr(datetime(T2,'convertfrom','juliandate'),'mmm-yyyy');

% Massen der Himmelskörper in kg
DiaST = 26000;    %Durchmesser Swift-Tuttle in m
rhoST = 2.5;      %Dichte Swift-Tuttle in kg/m^3
mST   = rhoST*4*pi*DiaST^3/3;  %Masse Swift-Tuttle in kg

mP = [2.2031868550e+13, 3.2485859200e+14, 3.9860043544e+14, ...
    4.2828375214e+13, 1.2668653190e+17,  3.7931206234e+16, ...
    5.7939512560e+15, 6.8350999700e+15,...  %7.2e13*G, 
    mST*G, ...
    1.3271244004e+20]/ G;
summP= sum(mP);                   % ungefähre Masse des Sonnensystems
epsE = EpsErde(T1);               % Ekliptikschiefe @ T1 


%% Keplerlösung und Reihenentwicklung
% Einlesen der Bahnparameter und Berechnung nach Keplerlösung 
% Achtung Rechnung hier in AE, AE/d 
Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T1,Aequi);
for indexP=1:PlanetsOnly
    PlanetKep(indexP) =PlanetPQR(T_vector, BaPa, BaPadot, indexP);
end

% Einlesen der Parameter und Berechnung der Erdbahn unter Berücksichtigung
% Mond und anderer Störungen (Reihenentwicklung)
% Achtung Rechnung hier in AE, AE/d 
%Exakte Koordinaten aus der Störungstheorie
t1Jh = Jd2JJht(T_vector);   %Umrechnung Julianische Jahrhunderte
SunData  = PerturbImport('SonnePos.csv');
SunPos   = SonneExakt(t1Jh,SunData,epsE,'AE');
PlanetKep(3).xyz  = -SunPos.xyz;

%%-------------------------------------------------------------------------
% AnfangsPositionen [km] und Geschwindigkeiten [km/s] der Himmelskörper 
% aus : https://ssd.jpl.nasa.gov/horizons/
vecr0 = AE* [
    [+3.5044731459e-01, -3.7407373056e-02, -3.6089929256e-02];%Merkur
    [-7.6445799950e-02, +7.1938215866e-01, +1.3916787662e-02];%Venus
    [-1.8324185857e-01, +9.7104012872e-01, +1.2685882047e-04];%Erde
    [-8.7535457155e-01, -1.2654778061e+00, -5.1567061263e-03];%Mars
    [+4.6495009432e+00, -1.7912164410e+00, -9.6589634761e-02];%Jupiter
    [+6.9514956851e+00, -7.0632684291e+00, -1.5395337276e-01];%Saturn
    [+1.4389044567e+01, +1.3482065167e+01, -1.3633986727e-01];%Uranus
    [+2.9624686966e+01, -4.0872817514e+00, -5.9856130334e-01];%Neptun
%     [-1.4276620837e+00, -8.4076159551e-01, +1.8789874266e-01];%9P/Temple 1
%     [-3.9399052555e-01, +7.1652579640e-01, +1.1583105373e-01];%2010TK7     
    [-3.0976406209e+01   +1.4326506949e+01  -2.1315103526e+01];%Swift-Tuttle
    [-8.5808349915e-03, +3.3470429582e-03, +1.7309053212e-04]];%Sonne
% in m
vecv0 = AE / Tag * [
    [-2.2707582290e-03, +2.9204389319e-02, +2.5953443972e-03];%Merkur
    [-2.0208176480e-02, -2.0266237759e-03, +1.1383547310e-03];%Venus
    [-1.7213898896e-02, -3.1295322660e-03, +3.5869599301e-07];%Erde
    [+1.2076503589e-02, -6.7024690766e-03, -4.3646377051e-04];%Mars
    [+2.6221777732e-03, +7.3957409953e-03, -8.9347514907e-05];%Jupiter
    [+3.6639554949e-03, +3.9017145704e-03, -2.1374264507e-04];%Saturn
    [-2.7180766241e-03, +2.6868888270e-03, +4.5192378845e-05];%Uranus
    [+4.0824560317e-04, +3.1283037388e-03, -7.3829594992e-05];%Neptun
%     [+1.0774018377e-02, -1.1795045602e-02, -2.6462859590e-03];%9P/Temple 1
%     [-1.6051208207e-02, -1.1225888345e-02, +6.5703843151e-03];%2010TK7
    [-1.4785875645e-03  +9.8592117066e-04  -4.8471485678e-04]; %Swift Tuttle
    [-3.3551917266e-06, -8.4435230812e-06, +1.4516419864e-07]];%Sonne
% in m/s

% Schwerpunktberechnung
vecrs = [0 0 0]; vecvs = [0 0 0];
for k=1:NK
    vecrs = vecrs + mP(k).*vecr0(k,:);
    vecvs = vecvs + mP(k).*vecv0(k,:);
end
vecrs = vecrs/summP;
vecvs = vecvs/summP;
for k=1:NK
    vecr0(k,:) = vecr0(k,:) - vecrs;
    vecv0(k,:) = vecv0(k,:) - vecvs;
end

%% Berechnung Trajektorien 

AB = [];
for k=1:NK
  AB = horzcat(AB,vecr0(k,:));
end
for k=1:NK
  AB = horzcat(AB,vecv0(k,:));
end
P1.G   = G;
P1.mP  = mP;
P1.NK  = NK;
% Solver ode45
opts = odeset('AbsTol',1.e-6,'RelTol',1.e-5);  %schnelle Übersichtsrechnung
% opts = odeset('AbsTol',1.e-10,'RelTol',1.e-11);%genaue Rechnung
disp('  Computation running !');
[tA,YA]=ode45(@(t,YA, P1)DGL_NBodySystem(t,YA,P1),t_vector,AB,opts,P1); 
disp('  Completed !');


%% Berechnung der heliozentrischen Koordinaten und Normierung auf AE
for achse =1:3
    SunCalc.xyz(:,achse)=YA(:,(NK-1)*3+achse);
end
for indexP=1:NK
    kS=(indexP-1)*3;
    for achse =1:3
      PlanetCalc(indexP).xyz(achse,:)=(YA(:,kS+achse)- ...
              SunCalc.xyz(:,achse))/AE;
    end
    PlanetCalc(indexP).r(:) = sqrt(PlanetCalc(indexP).xyz(1,:).^2 + ...
    PlanetCalc(indexP).xyz(2,:).^2 + PlanetCalc(indexP).xyz(3,:).^2);
end

dtplot= datetime(T1+tA/Tag,'convertfrom','juliandate');

Komet = PlanetCalc(NK-1);

%Berechnung Abstand Erde-Komet
rek=Komet.xyz-PlanetCalc(3).xyz;
rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));
[minrekbetr,ErdNaehe] = min(rekbetr);

%Berechnung Abstand Sonne-Komet und Perihel
rsk=Komet.xyz ;
rskbetr=sqrt(rsk(1,:).*rsk(1,:)+rsk(2,:).*rsk(2,:)+rsk(3,:).*rsk(3,:));
[minrskbetr,Perihel] = min(rskbetr);


%% Graphische Ausgabe  ----------------------------------------------------
header1 = 'Komet Swift-Tuttle';
header2 = strjoin([string(dt1str),' bis ',string(dt2str)]);
header2 = strjoin(['Swift-Tuttle: ',header2]);


% -------------------------------------------------------------------------
% Figure   Simulation Trajektorien 
% -------------------------------------------------------------------------

figure('name','Simulation der Bewegung der Planeten und des Kometen')
hold on
% Bahnen aus Keplerlösung';
for iPlot = 1:PlanetsOnly
   plot(PlanetKep(iPlot).xyz(1,:),...
         PlanetKep(iPlot).xyz(2,:),'linestyle',Style(1),...
         'Color',Colors(15,:),'LineWidth',1);
end
axis([-50 50 -50 50])
axis square
xlabel('x in AE')
ylabel('y in AE')
grid on
grid minor
hold on
%Bahnen aus Bewegungsgleichungen
for indexP=5:PlanetsOnly
    hp(indexP)=plot(PlanetCalc(indexP).xyz(1,:),PlanetCalc(indexP).xyz(2,:), ...
        'linestyle',Style(1),'color',Colors(indexP,:),'linewidth',1);
end
plot(Komet.xyz(1,:),Komet.xyz(2,:), ...
        'linestyle',Style(1),'color',Colors(9,:),'linewidth',1);
ht = title(strjoin({'Zeitraum : ',dt1str,'-',dt2str}));
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);
%Bahnen aus Bewegungsgleichungen
for k=1:366:NSteps
    for indexP=5:PlanetsOnly
    hs(indexP)=plot(PlanetCalc(indexP).xyz(1,k),PlanetCalc(indexP).xyz(2,k),'o', ...
         'markerfacecolor',Colors(indexP,:));
    end
    hk=plot(Komet.xyz(1,k),Komet.xyz(2,k), '+',...
        'color',Colors(9,:),'linewidth',2);
    Datumsausgabe = datestr(dtplot(k));
    htext = text(20,45, Datumsausgabe);
    pause(0.15);
    for indexP = 5:PlanetsOnly
        hs(indexP).Visible = 'off';
    end
    hk.Visible = 'off';
    htext.Visible = 'off';
end


% -------------------------------------------------------------------------
% Figure  1 Darstellung Trajektorien 
% -------------------------------------------------------------------------

figure('name','Bahnen der Planeten und des Kometen')
hold on
% Bahnen aus Keplerlösung';
for iPlot = 1:PlanetsOnly
   plot(PlanetKep(iPlot).xyz(1,:),...
         PlanetKep(iPlot).xyz(2,:),'linestyle',Style(1),...
         'Color',Colors(15,:),'LineWidth',1);
end
%Bahnen aus Bewegungsgleichungen
for indexP=1:PlanetsOnly
    plot(PlanetCalc(indexP).xyz(1,:),PlanetCalc(indexP).xyz(2,:), ...
        'linestyle',Style(1),'color',Colors(indexP,:),'linewidth',2);
end
plot(Komet.xyz(1,:),Komet.xyz(2,:), ...
        'linestyle',Style(1),'color',Colors(9,:),'linewidth',2);

%Anfangswerte Berechnung
for k=1:PlanetsOnly
    hp(k) =plot(vecr0(k,1)/AE, vecr0(k,2)/AE,'o','color',Colors(k,:),...
        'markerfacecolor',Colors(k,:),'markersize',4,'linewidth',2);
end
set(hp(5),'markersize',6)
hp(NK-1) =plot(vecr0(NK-1,1)/AE, vecr0(NK-1,2)/AE,'o','color',Colors(9,:),...
        'markerfacecolor',Colors(9,:),'markersize',4,'linewidth',2);
hp(NK)= plot(vecr0(NK,1)/AE, vecr0(NK,2)/AE,'o','color',Colors(10,:), ...
    'markersize',7,'markerfacecolor',Colors(10,:),'linewidth',3);
plot(0,0,'+k','markersize',10,'linewidth',1);  %Schwerpunkt
grid on
legend(hp, PNamen(1:end),'location', 'westoutside')
legend box off
axis([-50 50 -50 50])
axis equal
xlabel('x in AE')
ylabel('y in AE')
ht = title(strjoin({'Zeitraum : ',dt1str,'-',dt2str}));
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);


% -------------------------------------------------------------------------
% Figure 2 Darstellung Trajektorien im Zeitbereich um Perihel 1992 o. 2126 
% -------------------------------------------------------------------------
dt1P = datetime('2126-05-01 00:00:00');
dt2P = datetime('2126-09-30 00:00:00');
T1P  = juliandate(dt1P);
T2P  = juliandate(dt2P);

PNamen2    = PNamen(1:4);
PNamen2(5) = PNamen(9);
PNamen2(6) = PNamen(10);
PNamen2(7) = 'Erdnähe';
PNamen2(8) = 'Perihel';

dt1Pstr= datestr(datetime(T1P,'convertfrom','juliandate'),'dd.mm.');
dt2Pstr= datestr(datetime(T2P,'convertfrom','juliandate'),'dd.mm.');
plBeg = T1P -T1;
plEnd = T2P -T1;

figure('name','Bahnen der Planeten')
axis([-1 1.5 -1.5 1])
axis equal
hold on
msize = 6; %Markersize
% Bahnen aus Keplerlösung';
for iPlot = 1:4
   plot(PlanetKep(iPlot).xyz(1,:),...
         PlanetKep(iPlot).xyz(2,:),'linestyle',Style(1),...
         'Color',Colors(15,:),'LineWidth',1);
end
%Bahnen aus Bewegungsgleichungen nur bis Erde
for k=1:4 
    hp2(k) = plot(PlanetCalc(k).xyz(1,plBeg:plEnd),...
             PlanetCalc(k).xyz(2,plBeg:plEnd),'linestyle',Style(1),...
         'color',Colors(k,:),'linewidth',2);
end

%das Gleiche für Kometen
    hp2(5) = plot(Komet.xyz(1,plBeg:plEnd),...
         Komet.xyz(2,plBeg:plEnd),'linestyle',Style(1),...
         'color',Colors(9,:),'linewidth',2);

if T2P <= T2 
% Plot Zeitdaten
TimeStep = 30;
dtKstr = datestr(datetime(tA(plBeg:TimeStep:plEnd,1)/Tag+T1,'ConvertFrom',...
                'juliandate'),' dd.mm.');
dtErdNaehe = datestr(datetime(tA(ErdNaehe)/Tag+T1,'convertfrom','juliandate'),...
               ' dd.mm.yyyy ');
dtPerihel = datestr(datetime(tA(Perihel)/Tag+T1,'convertfrom','juliandate'),...
               ' dd.mm.yyyy ');
% Komet 
XLabelK = Komet.xyz(1,plBeg:TimeStep:plEnd);
YLabelK = Komet.xyz(2,plBeg:TimeStep:plEnd);
ZLabelK = Komet.xyz(3,plBeg:TimeStep:plEnd);
plot(XLabelK, YLabelK,'+', ...
     'color',Colors(9,:),'markerfacecolor',Colors(9,:),...
     'markersize',msize,'linewidth',2);
LabelPoints(XLabelK,YLabelK,dtKstr,'W',0.1,0.1,'FontSize',9,...
            'Color', Colors(9,:));
% Erde 
XLabel = PlanetCalc(3).xyz(1,plBeg:TimeStep:plEnd);
YLabel = PlanetCalc(3).xyz(2,plBeg:TimeStep:plEnd);
ZLabel = PlanetCalc(3).xyz(3,plBeg:TimeStep:plEnd);
plot(XLabel, YLabel,'+', ...
     'color',Colors(3,:),'markerfacecolor',Colors(3,:),...
     'markersize',msize,'linewidth',2);
LabelPoints(XLabel,YLabel,dtKstr,'E',0.1,0.1,'FontSize',9,...
            'Color', Colors(3,:));
%Sonne
hp2(6)= plot(0, 0,'o','color',Colors(10,:), ...
'markersize',msize+2,'markerfacecolor',Colors(10,:),'linewidth',2);
%Erdnähe
hp2(7) = line([PlanetCalc(3).xyz(1,ErdNaehe) ...
                  Komet.xyz(1,ErdNaehe)],...
                 [PlanetCalc(3).xyz(2,ErdNaehe) ...
                  Komet.xyz(2,ErdNaehe)], ...
                  'color', Colors(4,:), 'Linewidth',3);
LabelPoints(Komet.xyz(1,ErdNaehe),Komet.xyz(2,ErdNaehe), ...
            dtErdNaehe,'W',0.1,0.1,'FontSize',9,'Color', Colors(4,:));
hp2(8) = line([0 Komet.xyz(1,Perihel)],[0 Komet.xyz(2,Perihel)], ...
               'color', Colors(10,:), 'Linewidth',3);
LabelPoints(Komet.xyz(1,Perihel),Komet.xyz(2,Perihel), ...
            dtPerihel,'W',0.1,0.1,'FontSize',9,'Color', Colors(10,:));
grid on
legend(hp2(1:8), PNamen2(1:8), 'location', 'westoutside')
legend box off
axis equal
axis([-1 1.5 -1.5 1])
xlabel('x in AE')
ylabel('y in AE')
ht = title(strjoin({'Zeitraum : ',dt1Pstr,' bis ',dt2Pstr}));
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);
else
    disp('   ');
    disp('  Für Perihel 2126 bitte Integrationsreichweite vergrößern !');  
end

% -------------------------------------------------------------------------
% Figure 3 3D-Darstellung im Zeitbereich um Perihel 1992 o. 2126 
% -------------------------------------------------------------------------

figure('Name',header1);
xAx = 2;
plotrangeAll = plBeg:plEnd;

hp2(5)=plot3(Komet.xyz(1,plotrangeAll),Komet.xyz(2,plotrangeAll), ...
        Komet.xyz(3,plotrangeAll),'Color', ...
        Colors(9,:),'LineWidth',2);
hold on
axis equal
% Achse Sonne-Fruehlingspunkt 
SonneFP;
for iPlot= 1:4
  hp2(iPlot) = plot3(PlanetKep(iPlot).xyz(1,plotrangeAll), ...
      PlanetKep(iPlot).xyz(2,plotrangeAll), ...
      PlanetKep(iPlot).xyz(3,plotrangeAll),'Color', Colors(iPlot,:), ...
      'linewidth',2);
  plot3(PlanetKep(iPlot).xyz(1,:),PlanetKep(iPlot).xyz(2,:), ...
      PlanetKep(iPlot).xyz(3,:),'Color', Colors(iPlot,:), ...
      'linewidth',1);
end
hp2(6)= plot(0, 0,'o','color',Colors(10,:), ...
'markersize',msize+2,'markerfacecolor',Colors(10,:),'linewidth',2);

TimeStep = 30;
msize    = 10;
labelrange = plBeg:TimeStep:plEnd;
XLabel   = PlanetKep(3).xyz(1,labelrange);
YLabel   = PlanetKep(3).xyz(2,labelrange);
ZLabel   = PlanetKep(3).xyz(3,labelrange);
XLabelK  = Komet.xyz(1,labelrange);
YLabelK  = Komet.xyz(2,labelrange);
ZLabelK  = Komet.xyz(3,labelrange);

dtKstr = datestr(datetime(PlanetKep(3).Time(labelrange),'ConvertFrom',...
                'juliandate'),' dd.mm.');

%Zeitdaten
plot3(XLabel, YLabel, ZLabel,'+', ...
     'color',Colors(3,:),'markerfacecolor',Colors(9,:),...
     'markersize',msize,'linewidth',2);
text(XLabel, YLabel, ZLabel, dtKstr, 'FontSize',9,...
            'Color', Colors(3,:));

plot3(XLabelK, YLabelK, ZLabelK, '+', ...
     'color',Colors(9,:),'markerfacecolor',Colors(9,:),...
     'markersize',msize,'linewidth',2);
text(XLabelK,YLabelK, ZLabelK, dtKstr,'FontSize',9,...
            'Color', Colors(9,:));

hp2(7)=line([PlanetCalc(3).xyz(1,ErdNaehe) Komet.xyz(1,ErdNaehe)],...
     [PlanetCalc(3).xyz(2,ErdNaehe) Komet.xyz(2,ErdNaehe)],...
     [PlanetCalc(3).xyz(3,ErdNaehe) Komet.xyz(3,ErdNaehe)],...
     'Color', Colors(4,:),'Linewidth',3);
plot3(PlanetCalc(3).xyz(1,ErdNaehe), ...
    PlanetCalc(3).xyz(2,ErdNaehe),PlanetCalc(3).xyz(3,ErdNaehe), ...
    'o','Color', Colors(iPlot,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(4,:));
plot3(Komet.xyz(1,ErdNaehe), ...
    Komet.xyz(2,ErdNaehe), Komet.xyz(3,ErdNaehe), ...
    'o','Color', Colors(iPlot,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(4,:));
text(Komet.xyz(1,ErdNaehe), ...
     Komet.xyz(2,ErdNaehe), Komet.xyz(3,ErdNaehe),...
     dtErdNaehe, 'FontSize',9,'Color', Colors(4,:));

hp2(8)=line([0 Komet.xyz(1,Perihel)],...
     [0 Komet.xyz(2,Perihel)],...
     [0 Komet.xyz(3,Perihel)],...
     'Color', Colors(10,:),'Linewidth',3);
plot3(Komet.xyz(1,Perihel), ...
    Komet.xyz(2,Perihel), Komet.xyz(3,Perihel), ...
    'o','Color', Colors(10,:),'MarkerSize',5, ...
    'MarkerFaceColor',Colors(10,:));
text(Komet.xyz(1,Perihel), ...
     Komet.xyz(2,Perihel), Komet.xyz(3,Perihel),...
     dtPerihel, 'FontSize',9,'Color', Colors(10,:));


% Einstellungen der Achsen 
xlim([-xAx xAx]);
ylim([-xAx xAx]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));
shading flat
alpha 0.1
grid on;
axis equal;
legend(hp2, PNamen2, 'location', 'westoutside')
legend box off
axL = 1.6;
xlim([-axL axL]);
ylim([-axL axL]);
zlim([-axL axL]);
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');
ht = title(strjoin({'Zeitraum : ',dt1Pstr,'-',dt2Pstr}));
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);

%% Ausgabe Werte für Erdnähe 

dtErdNaehe = datestr(datetime(PlanetKep(3).Time(ErdNaehe),'convertfrom','juliandate'),...
               ' dd-mmm-yyyy HH:MM ');
dtPerihel  = datestr(datetime(PlanetKep(3).Time(Perihel),'convertfrom','juliandate'),...
               ' dd-mmm-yyyy HH:MM');

fprintf('\n')
fprintf('\n')
fprintf('\n Berechnung für 2126 über numerische Lösung der Bewegungsgleichungen')
fprintf('\n')
fprintf('\n Abstand Erdnähe : %12.8f AE ',minrekbetr)
fprintf('\n Abstand Erdnähe : %12.4e km ',minrekbetr*AEk)
fprintf('\n Datum Erdnähe   : %s', dtErdNaehe)
fprintf('\n')
fprintf('\n Abstand Perihel : %12.8f AE ',minrskbetr)
fprintf('\n Abstand Perihel : %12.4e km ',minrskbetr*AEk)
fprintf('\n Datum Perihel   : %s', dtPerihel)
fprintf('\n')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% DGL_Dystem

function dY = DGL_NBodySystem(~, Y, P1)
    G   = P1.G;
    mP  = P1.mP;
    NK  = P1.NK;
    r = zeros(NK,3);
    for k=1:NK
        kr = (k-1)*3;
        r(k,:) = [Y(1+kr) Y(2+kr) Y(3+kr)];
    end 
    a = zeros(NK,3);
    for iz=1:NK  
        for jz=1:iz
            if jz ~= iz 
                dr = r(jz,:) - r(iz,:);
                g = G*dr./ vecnorm(dr).^3;
                a(iz,:) = a(iz,:) + g * mP(jz);
                a(jz,:) = a(jz,:) - g * mP(iz);
            end
        end
    end
    dY = [];
    for k=1:NK*3
       dY = vertcat(dY, Y(NK*3+k));
    end 
    for k=1:NK
       dY = vertcat(dY, a(k,1));
       dY = vertcat(dY, a(k,2));
       dY = vertcat(dY, a(k,3));
    end
end


%% Ende Programm


