% -------------------------------------------------------------------------
% GravityAssist01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet verschieden AbhÃ¤ngigkeiten eines  
% Gravity-Assist-Manövers am Jupiter
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


% Parameter hier lles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % m
mP  = [0.33022; 4.8685; 5.9737; 0.64185; 1898.7; 568.51; ...
       86.849; 102.44]*1e24;  % in kg
aP  = [0.3870893; 0.72333199; 1.000000; 1.52366231; 5.20336301;...
       9.53707032; 19.191264; 30.068964]*AE;  %in m
aJ  = aP(5);    
MJ  = mP(5);
ME  = mP(3);
MS  = 1.989e30;

% Planetendaten ab hier alles in km,s,kg
vP  = sqrt(G*MS./aP)/1000; % Planetengeschwindigkeiten heliozentrisch 
GS  = MS*G*1e-9;   % Sonne
GE  = ME*G*1e-9;   % Erde
GJ  = MJ*G*1e-9;   % Jupiter
GP  = mP*G*1e-9;   % alle Planeten 

rP  = [200,500,1000,2000,5000]*1e3;  % Minimalabstand zum Planeten in km
v1  = linspace(5,15,5);              % "Einfluggeschwindigkeit" in km/s

gamma1d= wrapTo360(linspace(0,180,91));   % Winkel zw. vPVec und v1Vec
gamma1 = deg2rad(gamma1d);


%% Berechnung Streuwinkel und Geschwindigkeitszuwachs als Funktion von rmin
%  (alles in km /s)

% Ausgansgwerte
% oBdA setzen wir Vektor vPVec als [-vp,0]

vPVec  = [-vP(5); 0];  %Beispiel Jupiter
vPJ    = vecnorm(vPVec);
v1     = v1(3);
v1Vec  = [-v1*cos(gamma1); v1*sin(gamma1)];  % v1 = xxx km/s

% Berechnung Geschwindigkeit v1' = v1s im System Sigma'
v1sVec = v1Vec -vPVec;
vinf    = vecnorm(v1sVec);

% Berechnung Streuwinkel theta' = thetas und Delta vs als Funktion von rmin
vF = sqrt(2*GJ./rP); % "Fluchtgeschwindigkeit" von Minimaldistanz
Str1 = "r_{P} = ";
Str3 = "v_{1} = ";
Str4 = string(num2str(v1,'%4.1f km/s'));
Str5 = "v_{P} = ";
Str6 = string(num2str(vPJ,'%4.1f km/s'));

for k=1:length(rP)
    thetas(k,:)  = 2*asin(1./(1+2*vinf.^2/vF(k)^2));
    thetasd(k,:) = rad2deg(thetas(k,:));
    deltavs(k,:) = vinf./(1+2*vinf.^2/vF(k)^2)/vPJ;
    Str2 = string(num2str(rP(k),'%6.0e km'));
    lgdstr(k,:)  = strcat(Str1,Str2);
end 
lgdstr(length(rP)+1,:)  = strcat(Str3,Str4);
lgdstr(length(rP)+2,:)  = strcat(Str5,Str6);


%% Berechnung Geschwindigkeit v2 als Funktion von v1
%  (alles in km /s)

% Berechnung Geschwindigkeit v2' = v2s im System Sigma'
% Berechnung Geschwindigkeit v2 im System Sigma

for k=1:length(rP)
    costhe(k,:) = cos(thetas(k,:));
    sinthe(k,:) = sin(thetas(k,:));
    v2sVecX(k,:) = costhe(k,:).*v1sVec(1,:) - sinthe(k,:).*v1sVec(2,:);
    v2sVecY(k,:) = sinthe(k,:).*v1sVec(1,:) + costhe(k,:).*v1sVec(2,:);
    v2VecX(k,:)  = v2sVecX(k,:) + vPVec(1);
    v2VecY(k,:)  = v2sVecY(k,:) + vPVec(2);
    v2(k,:)      = sqrt(v2VecX(k,:).*v2VecX(k,:)+v2VecY(k,:).*v2VecY(k,:));
    gamma2d(k,:) = acosd((v2(k,:).^2+vPJ^2-vinf.^2)./(2*vPJ.*v2(k,:)));
end



%% Graphische Ausgabe


% Graphik  # 1
figure('name', 'Streuwinkel, Delta v')
subplot(1,2,1)
hold on
for k=1:length(rP)
    h(k) = plot(gamma1d, thetasd(k,:), ...,
           'color', Colors(k,:), 'linestyle', Style(1), 'linewidth', 1);
end
legend(h, lgdstr(1:length(rP),:),'location','northeast');
legend box off
grid on
xlabel('\gamma_1 in °'); ylabel('\vartheta´ in °')
xticks([0 30 60 90 120 150 180])
yticks([0 30 60 90 120 150 180])
xlim([0 180]); ylim([0 180]);
ttl=title('Streuwinkel');
text(45,170,lgdstr(end-1,:),'FontSize',12,'FontWeight','normal');
text(45,160,lgdstr(end,:),'FontSize',12,'FontWeight','normal');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');

subplot(1,2,2)
hold on
for k=1:length(rP)
    h(k) = plot(gamma1d, v2(k,:)-v1, ...,
           'color', Colors(k,:), 'linestyle', Style(1), 'linewidth', 1);
end
legend(h, lgdstr(1:length(rP),:),'location','northeast');
legend box off
grid on
xlabel('\gamma_1 in °'); ylabel('v_2 - v_1  in km/s')
xlim([0 180]); xticks([0 30 60 90 120 150 180])
text(45,170,lgdstr(end-1,:),'FontSize',12,'FontWeight','normal');
text(45,160,lgdstr(end,:),'FontSize',12,'FontWeight','normal');
ttl=title('Geschwindigkeitsgewinn');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');



figure('name', 'Delta v, gamma2 - gamma1')
subplot(1,2,1)
hold on
for k=1:length(rP)
   h2(k) = plot(gamma1d, deltavs(k,:),'color', Colors(k,:),...
        'linestyle', Style(3),'linewidth', 1);
end
grid on
ylabel('\Delta v / v_P  '); xlabel('\gamma_1 in °')
legend(h2, lgdstr(1:length(rP),:),'location','northwest');
legend box off
xlim([0 180]); xticks([0 30 60 90 120 150 180])
ttl=title('\Delta v / v_P im Planetensystem');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');


subplot(1,2,2)
hold on
for k=1:length(rP)
   h2(k) = plot(gamma1d,gamma2d(k,:)-gamma1d,'color', Colors(k,:),...
        'linestyle', Style(3),'linewidth', 1);
end
grid on
ylabel('\gamma_2 - \gamma_1'); xlabel('\gamma_1 in °')
legend(h2, lgdstr(1:length(rP),:),'location','northwest');
legend box off
xlim([0 180]); xticks([0 30 60 90 120 150 180])
ttl=title('\gamma_2 - \gamma_1');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');

%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

function yout=R2D_z(phi)
% Eingabe:  
%   phi     = Rotationswinkel im BogenmaÃ [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_z 
   yout=[cos(phi) sin(phi);...
        -sin(phi) cos(phi)];
end

%Ende Funktionen
% -------------------------------------------------------------------------

