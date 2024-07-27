%%-------------------------------------------------------------------------
% SonnensystemSimulationVergleichJPL.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Himmelsmechanik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Simulation der Bewegung der Planeten des Sonnensystems
% um den Schwerpunkts des Sonnensystem (Solar System Barycenter - SSBC)
% und Vergleich der numerischen Ergebnisse mit JPL-Daten, sowie den
% Rechnungen in O. Natt, Physik mit Python, Springer Verlag 2022
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
          "Neptun" "9P/Tempel 1" "Swift-Tuttle" "Sonne"];
NK = length(PNamen);


% Parameter, Rechnungen erfolgen alle in m/s
AE   = 1.495978707e11;            % AE in m
G    = 6.6743e-11;                % Gravitationskonstante m^3 / kg / s^2
Tag  = 24 * 60 * 60;              % in s
RSun = 696342e3; rsun = RSun/AE;  % Sonnenradius in m bzw AE
RE   = 6378000;  rE   = RSun/AE;  % Erdradius in m bzw AE

%% Einladen der JPL Daten Erde und Mars für 2022

% Daten in den csv-Files sind in km und km/s
% Umrechnung in m und m/s
JPLD = JPLVectorTable('JPL_Earth_SSB_2022.csv');
NSteps = length(JPLD.T);  % In Schritten von 1h
for k= 1:NSteps
    JPLE.r(k) = sqrt((JPLD.x(k))^2+(JPLD.y(k))^2+ (JPLD.z(k))^2)*1000;
end
JPLD = JPLVectorTable('JPL_Mars_SSB_2022.csv');
for k= 1:NSteps
    JPLM.r(k) = sqrt((JPLD.x(k))^2+(JPLD.y(k))^2+ (JPLD.z(k))^2)*1000;
end

% Massen der Himmelskörper in kg
DiaST = 26000;    %Durchmesser Swift-Tuttle in m
rhoST = 1.5;      %Dichte Swift-Tuttle in kg/m^3
mST   = rhoST*4*pi*DiaST^3/3;  %Masse Swift-Tuttle in kg

mP = [2.2031868550e+13, 3.2485859200e+14, 3.9860043544e+14, ...
    4.2828375214e+13, 1.2668653190e+17,  3.7931206234e+16, ...
    5.7939512560e+15, 6.8350999700e+15,  7.2e13*G, mST*G, ...
    1.3271244004e+20]/ G;
summP= sum(mP);

% Simulationszeit und Zeitschrittweite [in Tagen].
dt1str= "2022-Jan-01";
dt1   = datetime(dt1str);    
T1 = juliandate(dt1); % Julianisches Datum
dt1str= datestr(datetime(T1,'convertfrom','juliandate'),'mmm-yyyy');

NYears = 1;
dt    = 3600;              %Schrittweite in s
tend  = NSteps*dt;         %Ende der Simulation in s
Tend  = T1+(tend-dt)/Tag;  %Ende der Simulation in Julian. Datum
t_vector=0:dt:tend-dt;
dt2str = datestr(datetime(Tend,'convertfrom','juliandate'),'mmm-yyyy');

%%------------------------------------------------------------------------
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
    [-1.4276620837e+00, -8.4076159551e-01, +1.8789874266e-01];%9P/Temple 1
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
    [+1.0774018377e-02, -1.1795045602e-02, -2.6462859590e-03];%9P/Temple 1
    %[-1.6051208207e-02, -1.1225888345e-02, +6.5703843151e-03];%2010TK7
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

% Vergleichsdaten O. Natt
filename = 'ephemeriden.mat';
A= importdata(filename);
B=A.r(4,:);
C=A.r(5,:);
for k=1:NSteps
  kS = (k-1)*3+1;
  for achse =1:3
  E.xyz(k,achse) = B(kS+achse-1);
  end
  E.t(k) = A.t(k);
end
E.r = vecnorm(E.xyz');
for k=1:NSteps
  kS = (k-1)*3+1;
  for achse =1:3
  M.xyz(k,achse) = C(kS+achse-1);
  end
end
M.r = vecnorm(M.xyz');

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
opts = odeset('AbsTol',1.e-13,'RelTol',1.e-12);
disp('  Computation running !');
[tA,YA]=ode45(@(t,YA, P1)DGL_NBodySystem(t,YA,P1),t_vector,AB,opts,P1); 
disp('  Completed !');

% Berechnung Schwerpunktlage zum SSBC
RS = zeros(length(tA),3);
for indexP=1:NK
    kS = (indexP-1)*3;
    for pos=1:length(tA)
        for achse =1:3
            Planet(indexP).xyz(pos,achse)=YA(pos,kS+achse);
            RS(pos,achse) = RS(pos,achse) +...
                  mP(indexP).*Planet(indexP).xyz(pos,achse);
        end
        Planet(indexP).r(pos) = vecnorm(Planet(indexP).xyz(pos,:));
    end
end
RS =RS/summP;

% Ausgabe Werte Erde in Datenfile
fileID = fopen('./Data/EarthCalc2022.txt','w');
dtplot= datetime(T1+tA/Tag,'convertfrom','juliandate');
fprintf(fileID, '%s \t %s \t %s \t %s \n', 't', 'x', 'y', 'z');

for k=1:length(tA)
    fprintf(fileID, '%s %12.10e %12.10e %12.10e \n', ...
           datestr(dtplot(k),'dd-mmm-yyyy hh:MM'), ...
            Planet(3).xyz(k,1),Planet(3).xyz(k,2),Planet(3).xyz(k,3));
end
fclose(fileID);


%% Graphische Ausgabe Vergleich Berechnungen mit Kepler und JPL

for indexP=3:4
% indexP = 3 für Erde
% indexP = 4 für Mars
% Differenzen der Berechnungen (volle Integration, Kepler, JPL)
    if indexP == 3
           Diff1 = (Planet(indexP).r(:) - JPLE.r(:))/RE;
           Diff2 = Planet(indexP).r(:)  - E.r(:);
    else
           Diff1= (Planet(indexP).r(:)  - JPLM.r(:))/RE;
           Diff2 = Planet(indexP).r(:)  - M.r(:);
    end
    figure('Name',strcat('Vergleich - ',PNamen(indexP)))
    yyaxis left
    h(1)=plot(dtplot,Diff1,'Color',Colors(3,:),'LineWidth',2);
    ylabel('Delta r in R_E')
    if indexP == 3 
        ylim([-2 +2]);
    else 
        ylim([-0.2 0.2]);
    end
    yyaxis right
    hold on
    h(2)=plot(dtplot,Diff2,'Color',Colors(4,:),'LineWidth',2);
    ylabel('Delta r in m')
    if indexP == 3 
        ylim([-2 +2]);
    else 
        ylim([-1e-1 1e-1])
    end
    xlabel('Zeit')
    grid on
    hl=legend(h, 'Abweichung zu JPL-Daten',...
                 'Abweichung zur Python-Rechnung [O. Natt]', ...
                 'location', 'north');
    legend box off
    ht = title(strcat("Vergleich für ",PNamen(indexP)));
    set(ht,'FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',14)
end

%-------------------------------------------------------------------
% Darstellung Schwerpunkt-Bewegung

figure('Name','Schwerpunkt-Bewegung')
plot(dtplot(:,:),RS(:,1));
hold on
plot(dtplot(:,:),RS(:,2));
grid on
xlabel('Zeit')
ylabel('x_S,y_S in m')
legend('x','y', 'location', 'best');
legend box off
ht = title('Bewegung des Schwerpunktes');
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14)

%-------------------------------------------------------------------
% Darstellung Trajektorie Sonne

figure('name','Bahnen des SSBC')
hold on
%Bahnen aus Bewegungsgleichungen
indexP = NK;
kS = (indexP-1)*3;
drawcircle('Center',[-YA(1,kS+1)/AE,-YA(1,kS+2)/AE], ...
           'Radius',rsun,'Color',Colors(10,:));
hp(1)  = viscircles([-YA(1,kS+1)/AE,-YA(1,kS+2)/AE],rsun, ...
                    'Color',Colors(10,:));
viscircles([-YA(end,kS+1)/AE,-YA(end,kS+2)/AE],rsun,'Color',Colors(10,:));
hp(2)  = plot(0,0,'+k','markersize',10,'linewidth',2);  %Schwerpunkt
hp(3)  = plot(-YA(1:end,kS+1)/AE,-YA(1:end,kS+2)/AE,...
         'linestyle',Style(1),...
         'color','k','linewidth',2);
hp(4)  = plot(-YA(1:end,1)/AE,-YA(1:end,2)/AE,...
         'linestyle',Style(3),...
         'color',Colors(1,:),'linewidth',1);
grid on
hl=legend(hp, 'Größe der Sonne', 'Solar System Barycenter', ...
              'Bewegung Sonne', 'Bewegung Merkur','location', 'south');
legend box off
set(hl,'FontSize',14,'FontWeight','normal');
axis equal
axis([-0.02 0.02 -0.02 0.02])
ylabel('x in AE')
xlabel('y in AE')
ht = title(strjoin({'Zeitraum : ',dt1str,'-',dt2str}));
set(ht,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%%%%%%%%%%%%%
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


