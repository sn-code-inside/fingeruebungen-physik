% -------------------------------------------------------------------------
% GravityAssist02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet verschieden AbhÃ¤ngigkeiten beim 
% Gravity-Assist-Manöver.
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Parameter hier alles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % m
mP  = [0.33022; 4.8685; 5.9737; 0.64185; 1898.7; 568.51; ...
       86.849; 102.44]*1e24;  % in kg
aP  = [0.3870893; 0.72333199; 1.000000; 1.52366231; 5.20336301;...
       9.53707032; 19.191264; 30.068964]*AE;  %in m
NamesP = ["Merkur ","Venus  ","Erde   ","Mars   ","Jupiter",...
             "Saturn ","Uranus ","Neptun"];
aJ  = aP(5);    
MJ  = mP(5);
ME  = mP(3);
MS  = 1.989e30;

% Planetendaten ab hier alles in km,s,kg
vP  = sqrt(G*MS./aP)/1000; % Planetengeschwindigkeiten heliozentrisch 
GS  = MS*G*1e-9;   % Sonne
GP  = mP*G*1e-9;   % alle Planeten 

RP   = [2439,6051,6378,3394,71400,60000,25650,24780];%Planetenradius in km
RSOI = (aP.*(mP/MS).^(2/5))/1000;  %Radius SOI in km


gamma1d = [20,20,20,20,60,60,60,60];                        % Winkel zw. vPVec und v1Vec
gamma1  = deg2rad(gamma1d);
v1     = 12.7;
v1Vec  = [-v1.*cos(gamma1); v1.*sin(gamma1)];  % v1 = xxx km/s


%% Berechnung Streuwinkel und Geschwindigkeitszuwachs als Funktion von rmin
%  (alles in km /s)

% Ausgansgwerte
% oBdA setzen wir Vektor vPVec als [-vp,0]
zerVec = zeros(1,length(vP));
vPVec  = [-vP(:), zerVec' ]';  %Beispiel Jupiter
vPBet  = vecnorm(vPVec);

% Berechnung Geschwindigkeit v1' = v1s im System Sigma'
v1sVec = v1Vec -vPVec;
vinf   = vecnorm(v1sVec);

% Berechnung Streuwinkel theta' = thetas und Delta vs als Funktion von rmin

for k = 1:length(aP) %Schleife über Planeten (k)
    rP(k,:)=linspace(RP(k)+300,RSOI(k)/10,1000);
    vF(k,:) = sqrt(2*GP(k)./rP(k,:)); % "Fluchtgeschwindigkeit" von Minimaldistanz
end
Str3 = "v_{1} = ";
Str4 = string(num2str(v1,'%4.1f km/s'));
Str5 = "\gamma_{1} = ";
Str6 = string(num2str(gamma1d(1),'%4.1f °'));
Str7 = "\gamma_{1} = ";
Str8 = string(num2str(gamma1d(5),'%4.1f °'));

for k=1:4
   lgdstr1(k,:) = strcat(" \vartheta - ",NamesP(k)); 
   lgdstr1(k+4,:) = strcat(" \Delta v - ",NamesP(k));
   lgdstr2(k,:) = strcat(" \vartheta - ",NamesP(k+4)); 
   lgdstr2(k+4,:) = strcat(" \Delta v - ",NamesP(k+4));
end 

for k=1:length(aP) %Schleife über Planeten (k)
    thetas(k,:)  = 2*asin(1./(1+2*vinf(k).^2./vF(k,:).^2));
    thetasd(k,:) = rad2deg(thetas(k,:));
    deltavs(k,:) = vinf(k)./(1+2*vinf(k).^2./vF(k,:).^2)./vPBet(k);
end 
lgdstr(length(rP)+1,:)  = strcat(Str3,Str4);
lgdstr(length(rP)+2,:)  = strcat(Str5,Str6);
lgdstr(length(rP)+3,:)  = strcat(Str7,Str8);


%% Berechnung Geschwindigkeit v2 als Funktion von v1
%  (alles in km /s)

% Berechnung Geschwindigkeit v2' = v2s im System Sigma'
% Berechnung Geschwindigkeit v2 im System Sigma

for k=1:length(aP)
    costhe(k,:) = cos(thetas(k,:));
    sinthe(k,:) = sin(thetas(k,:));
    v2sVecX(k,:) = costhe(k,:).*v1sVec(1,k) - sinthe(k,:).*v1sVec(2,k);
    v2sVecY(k,:) = sinthe(k,:).*v1sVec(1,k) + costhe(k,:).*v1sVec(2,k);
    v2VecX(k,:)  = v2sVecX(k,:) + vPVec(1,k);
    v2VecY(k,:)  = v2sVecY(k,:) + vPVec(2,k);
    v2(k,:)      = sqrt(v2VecX(k,:).*v2VecX(k,:)+v2VecY(k,:).*v2VecY(k,:));
    gamma2d(k,:) = acosd((v2(k,:).^2+vP(k).^2-vinf(k).^2)./(2*vP(k).*v2(k,:)));
end



%% Graphische Ausgabe

% Graphik  # 1
figure('name', 'Streuwinkel, Delta v')
subplot(1,2,1)
yyaxis left
hold on
for k=1:4
    h(k) = plot(rP(k,:)/RSOI(k), thetasd(k,:), ...,
           'color', Colors(k,:), 'linestyle', Style(1), 'linewidth', 1);
end
xlabel('r_P / R_{SOI}'); ylabel('\vartheta´ in °')
text(0.08,12.0,lgdstr(end-2,:),'FontSize',12,'FontWeight','normal');
text(0.08,11.0,lgdstr(end-1,:),'FontSize',12,'FontWeight','normal');
yyaxis right
hold on
for k=1:4
    h2(k) = plot(rP(k,:)/RSOI(k), (v2(k,:)-v1), ...,
           'color', Colors(k,:), 'linestyle', Style(4), 'linewidth', 2);
end
xlabel('r_P / R_{SOI}'); ylabel('\Delta v = v_2 - v_1  in km/s')
grid on
legend([h,h2], lgdstr1(1:8,:) ,'location','northeast');
legend box off
ttl=title('Streuwinkel, \Delta v');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');

% Graphik  # 2
subplot(1,2,2)
yyaxis left
hold on
for k=5:8
    h(k-4) = plot(rP(k,:)/RSOI(k), thetasd(k,:), ...,
           'color', Colors(k,:), 'linestyle', Style(1), 'linewidth', 1);
end
xlabel('r_P / R_{SOI}'); ylabel('\vartheta´ in °')
text(0.08,95.0,lgdstr(end-2,:),'FontSize',12,'FontWeight','normal');
text(0.08,85.0,lgdstr(end,:),'FontSize',12,'FontWeight','normal');
yyaxis right
hold on
for k=5:8
    h2(k-4) = plot(rP(k,:)/RSOI(k), v2(k,:)-v1, ...,
           'color', Colors(k,:), 'linestyle', Style(4), 'linewidth', 2);
end
xlabel('r_P / R_{SOI}'); ylabel('\Delta v = v_2 - v_1  in km/s')
grid on
legend([h,h2], lgdstr2(1:8,:) ,'location','northeast');
legend box off
ttl=title('Streuwinkel, \Delta v');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');

%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

function yout=R2D_z(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_z 
   yout=[cos(phi) sin(phi);...
        -sin(phi) cos(phi)];
end

%Ende Funktionen
% -------------------------------------------------------------------------

