% -------------------------------------------------------------------------
% GravityAssistAnalytisch.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet verschieden AbhÃƒÂ¤ngigkeiten beim
% Gravity-Assist-ManÃ¶ver mit analytischen Formeln.
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


%% Parameter hier alles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % m
MS  = 1.989e30;

% Planetendaten ab hier alles in km,s,kg
mP  = [0.33022; 4.8685; 5.9737; 0.64185; 1898.7; 568.51; ...
       86.849; 102.44]*1e24;  % in kg
aP  = [0.3870893; 0.72333199; 1.000000; 1.52366231; 5.20336301;...
       9.53707032; 19.191264; 30.068964]*AE;  %in m
NamesP = ["Merkur ","Venus  ","Erde   ","Mars   ","Jupiter",...
             "Saturn ","Uranus ","Neptun "];
vP  = sqrt(G*MS./aP)/1000; % Planetengeschwindigkeiten heliozentrisch 
muS = MS*G*1e-9;   % Sonne in km^3/s^2
muP = mP*G*1e-9;   % alle Planeten 
RP   = [2439,6051,6378,3394,71400,60000,25650,24780];%Planetenradien in km
RSOI = (aP.*(mP/MS).^(2/5))/1000;  %Radien der SOI in km

% Beispiel Jupiter
aJ  = aP(5);    
MJ  = mP(5);
vJ  = vP(5);
RJ  = RP(5);
muJ = muP(5);

%% Berechnung v2 und Streuwinkel - Beispiel Jupiter

% Input
% Variation v1

gamma1d = linspace(10,170,161);   % Winkel zw. vPVec und v1Vec
gamma1  = deg2rad(gamma1d);
v1     = [11.0,13.0,15.0];  % v1 = xxx km/s
rPer   = 5*RJ;
for k=1:length(v1)
   v1S(k,:)  = sqrt(v1(k).^2+vJ^2-2*v1(k)*vJ.*cos(gamma1(:)));
   C1(k,:)   = v1(k)^2-vJ^2-v1S(k,:).^2;
   C2(k,:)   = 1+rPer*v1S(k,:).^2/muJ;
   v2q(k,:)  = v1(k)^2-2*C1(k,:)./C2(k,:).^2+(4*v1S(k,:)*vJ./C2(k,:)).*...
               sqrt(1-1./C2(k,:).^2).*...
               sqrt(1-(C1(k,:)/2./v1S(k,:)/vJ).^2);
   v2(k,:)   = sqrt(v2q(k,:));
end

figure(1)
hold on
for k=1:length(v1)
  h(k) = plot(gamma1d, v2(k,:));  
  lgdstr(k,:) = strjoin(["v_1 =",num2str(v1(k),'%4.1f km/s')]);
 end
grid on
ttl = title('Geschwindigkeit v_2 über \gamma_1' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('\gamma_1 in °','FontSize',13); 
ylabel('v_2 in km/s','FontSize',14);
lgd=legend(h,lgdstr,...
                 'location','south');
legend box off
set(lgd,'FontSize',14)
set(gca,'FontSize',14)


%% Berechnung v2 und Streuwinkel - Beispiel Jupiter

% Input
% Variation rPer

gamma1d = 60;   % Winkel zw. vPVec und v1Vec
gamma1  = deg2rad(gamma1d);
v1     = [11.0,13.0,15.0];  % v1 = xxx km/s
rPer   = linspace(3,100,98)*RJ;
for k=1:length(v1)
   v1S(k)   = sqrt(v1(k).^2+vJ^2-2*v1(k)*vJ.*cos(gamma1));
   C1a(k)   = v1(k)^2-vJ^2-v1S(k).^2;
   C2a(k,:) = 1+rPer(:).*v1S(k).^2/muJ;
   v2qa(k,:)= v1(k)^2-2*C1a(k)./C2a(k,:).^2+(4*v1S(k)*vJ./C2a(k,:)).*...
               sqrt(1-1./C2a(k,:).^2).*...
               sqrt(1-(C1a(k)/2./v1S(k)/vJ).^2);
   v2a(k,:) = sqrt(v2qa(k,:));
end


figure(2)
hold on
for k=1:length(v1)
  h2(k)=plot(rPer(:)/RJ, v2a(k,:));  
end
grid on
ttl = title('Geschwindigkeit v_2 über r_{per}' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('r_{per} in R_J','FontSize',13); 
ylabel('v_2 in km/s','FontSize',14);
lgd=legend(h2,lgdstr,...
                 'location','northeast');
legend box off
set(lgd,'FontSize',14)
set(gca,'FontSize',14)

%% Berechnung Energiegewinn bei rPer = 5*RP

rPeriz = 5*RJ;
beta1  = linspace(0,180,181);   % Winkel beta1
v1S    = linspace(5,25,181);    % in km/s
x      = v1S*sqrt(rPeriz/muJ);

for k=1:length(beta1)
  DelH(k,:)  = vJ*sqrt(muJ/rPeriz)*2.*x(:)./(1+x(:).^2).^2.*....
               (x(:).*sqrt(2+x(:).^2).*sind(beta1(k))-cosd(beta1(k)));
  DelH(k,:)  = DelH(k,:)/vJ/sqrt(muJ/rPeriz);
end

figure(3)
% contour(beta1(:),x(:), DelH(:,:),linspace(-1,1,20)); 
% hold on
% line([0 180],[1 1]);
% line([120 120],[0.2 1.4]);
surf(beta1(:),x(:), DelH(:,:)); 
xlim([0 180]); ylim([0.3 1.3]);
hold on
grid  on
ttl = title('Energiegewinn \Delta H' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('\beta_1 °','FontSize',13); 
ylabel('v_1 normiert auf v_F','FontSize',14);
set(gca,'FontSize',14)

gamma1dE = linspace(10,170,181);   % Winkel zw. vPVec und v1Vec
v1E      = linspace(10,25,181);    % in km/s

% Alternative Berechnung

for k=1:length(v1E)
   v1SE(k,:)  = sqrt(v1E(k).^2+vJ^2-2*v1E(k)*vJ.*cosd(gamma1dE(:)));
   C1E(k,:)   = v1E(k)^2-vJ^2-v1SE(k,:).^2;
   C2E(k,:)   = 1+rPeriz*v1SE(k,:).^2/muJ;
   v2qE(k,:)  = v1E(k)^2-2*C1E(k,:)./C2E(k,:).^2+(4*v1SE(k,:)*vJ./C2E(k,:)).*...
               sqrt(1-1./C2E(k,:).^2).*...
               sqrt(1-(C1E(k,:)/2./v1SE(k,:)/vJ).^2);
   DelHE(k,:) = 0.5*(v2qE(k,:)-v1E(k)^2)/vJ/sqrt(muJ/rPeriz);
end

figure(4)
contour(gamma1dE(:),v1E(:), DelHE(:,:)); 
% surf(gamma1dE(:),v1E(:), DelHE(:,:),linspace(-1,1,20)); 
hold on
v1max = sqrt(muJ/rPeriz);
line([10 170],[v1max v1max]);
line([75 75] ,[10 25]);
xlim([10 170]); ylim([10 25]);
grid  on
ttl = title('Energiegewinn \Delta H' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('\gamma_1 °','FontSize',13); 
ylabel('v_1 in km/s','FontSize',14);
set(gca,'FontSize',14)

%% Maximaler Energieübertrag Delta H . optimales  v1' und max theta'

%optimales v1' (durch v1 und gamma1 bestimmbar)
for k= 1:length(muP)
    v1Sopt(k) = sqrt(muP(k)./RP(k));
end

%max Delta H (durch v1 und gamma1 bestimmbar)
for k= 1:length(muP)
    DelHmax(k) = v1Sopt(k)*vP(k)*1e6;
end

%max theta'
v1Snormal = 12.7;
for k= 1:length(muP)
    vesc(k)     = sqrt(2*muP(k)/RP(k));
    thetamax(k) = 2*asind(1/(1+2*v1Snormal^2/vesc(k)/vesc(k)));
end


fprintf('\n')
fprintf('Parameter Gravity-Assist-Manöver Planeten unseres Sonnensystems:')
fprintf('\n') 
fprintf('|  Planet   |  RP (km)  |  v1opt (km/s) | Delta H (J/kg) | theta max ° |\n');

for k=1:length(muP)
fprintf('|  %s  |  %7.0f  |    %6.2f     |    %7.2e    |   %7.2f   |\n',...
                NamesP(k), RP(k), v1Sopt(k) , DelHmax(k), thetamax(k));
end
fprintf('\n')


%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

%Ende Funktionen
% -------------------------------------------------------------------------

