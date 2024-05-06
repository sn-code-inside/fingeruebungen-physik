% -------------------------------------------------------------------------
% LambertMarsVenus.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Parameter der Ellipsenbahn für Flüge zu Mars und
% Venus für das Lambert-Problem bei vorgegebenen Positionsvektoren 1 und 2 
% in AE, sowie der Abflug und Ankunftszeit. 
% -------------------------------------------------------------------------
%
% Berechnung Lambert-Problem
%
% Input:                                         Einheit
%   vecr1       - Positionsvektor 1                 AE
%   vecr2       - Positionsvektor 2                 AE
%   dm          - Richtung Bewegung                 'pro','retro'
%   TOF         - Flugzeit                          d
%
% OutPut (u.a.)
%   a           - Große Halbachse                   AE
%   p           - semi-latus rectum                 AE
%   exz         - Exzentrizität                  
%   ups0        - Periapsis                 
%   V1          - Geschwindigkeitsvektor 1         km/s
%   V2          - Geschwindigkeitsvektor 2         km/s
%
%
% Programm benutzt Routine LambertSolver3.m, die auf Basis der von
% Meysam Mahooti entwickelten Routine LambertBATTIN.m  erweitert wurde.
% Für die Routine LambertBATTIN.m liegt das 
% Copyright (c) bei Meysam Mahooti.
% (Siehe  LambertSolver3.m !)
% 
% -------------------------------------------------------------------------

clc
clear 
close all 
format long g

addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Randwertproblem 

AEk = 149597870.700;   % Astronomische Einheit in km
ds  = 86400;           % Tag in s

% Gravitationsparameter
muE  = 398600.4415;    % G*ME in km^3/s^2 für Erde
muS  = 2.95479E-04;    % AE^3/Tage^2 für Sonne

% Festlegung Zentralköper
ZKI = "Heimflug von ISS";   % ZK Erde,  Einheiten km,s
ZKV = "Flug zur Venus  ";   % ZK Sonne, Einheiten AE,s
ZKM = "Flug zum Mars   ";   % ZK Sonne, Einheiten AE,s
titlestr = [ZKV, ZKM, ZKI];
dm(1) = "pro";
dm(2) = "retro";
%Farben für Venus Darstellung
Co(1,3) = 10; Co(1,2) = 2; Co(1,1) = 3;
%Farben für Mars Darstellung
Co(2,3) = 10; Co(2,2) = 4; Co(2,1) = 3;
%Farben für Erdorbit Darstellung
Co(3,3) = 3; Co(3,2) = 3; Co(3,1) = 1;

% Abfrage der Eingabedaten in der Kommandozeile

fprintf('\nBeispiel-Probleme ');
Beisp = input ('\n[V] Venus\n[M] Mars\n[I] ISS\n[E] Exit .........: ','s');
Beisp = upper(Beisp);

%% Parameter Beispiel-Rechnungen
% Beispiel 1 Venusflug 
switch Beisp 
  case 'V'    
    ZK     = ZKV;
    mu     = muS;  
    ZKNr   = 1;
    r1     = 1.006;   % in AE
    r2     = 0.727;   % in AE
    gamma  = deg2rad(273.9) ;
    TOF    = 158;     % in d
    chi1   = deg2rad(-9.1);   %Anfangswinkel P1 in °
    chi2   = deg2rad(264.8);  %Endwinkel P2 in °
  case 'M' 
    ZK     = ZKM;
    mu     = muS;  
    ZKNr   = 2;
    r1     = 1.006;   % in AE
    r2     = 1.525;   % in AE
    gamma  = deg2rad(75) ;
    TOF    = 115;     % in d
    chi1   = deg2rad(30);     %Anfangswinkel P1 in °
    chi2   = chi1+gamma;      %Endwinkel P2 in °
   case 'I'    
    ZK     = ZKI;
    mu     = muE;  
    ZKNr   = 3;
    r1     = 6800;   % in km
    r2     = 6400;   % in km
    gamma  = deg2rad(75) ;
    TOF    = 3000;     % in s
    chi1   = pi-gamma/2;  %Anfangswinkel P2 in °
    chi2   = pi+gamma/2;  %Endwinkel P2 in °
end     
       
     
      
%% Lösung Lambert-Problem mit LambertSolver3.m


if Beisp =='V' ||  Beisp =='M'||  Beisp =='I' 
    vecr1  = r1*[cos(chi1),sin(chi1),0];
    vecr2  = r2*[cos(chi2),sin(chi2),0];
    V10  = sqrt(mu/r1)*[-sin(chi1),cos(chi1),0];
    if ZKNr < 3 
        V10=V10*AEk/ds;
    end
    
    for k=1:2
        % Solver gibt V1, V2 
        [V1, V2, a] = LambertSolver3(vecr1, vecr2, dm(k), TOF, mu);
        RLvec = (1/r1-1/a)*vecr1-(1/mu)*(dot(vecr1,V1))*V1; % Laplace-Vektor
        exz   = vecnorm(RLvec);
        p     = a*(1-exz^2);
        ups0 = atan2(RLvec(2),RLvec(1));                    
        cosE1  = (a-r1)/a/exz;
        ups1   = -acos(a*(cosE1-exz)/r1)+ups0;
        % Beträge und Umrechnung der Geschwindigkeiten in km/s
        if ZKNr < 3
             V1 = V1*AEk/ds;
             V2 = V2*AEk/ds;
        end
        v1 = vecnorm(V1);
        v2 = vecnorm(V2);
        % Winkel zwischen Ursprunsggeschwindigkeit und Post-Manöver-Geschw.
        % am Punkt P1
        delta = acos(dot(V1,V10)/vecnorm(V10)/v1);               
        % Ausdruck
        if delta < pi/2 
        % Wir betrachten nur die energetisch günstigere Variante, 
        % d.h. Winkel zwischen den Vektoren der Geschwindigkeit < +-90°,
        % d.h. nur die prograde Lösung.
            ttlprint = PrintOut(ZK, ZKNr, dm(k), r1, r2,...
                                TOF, gamma, ups0, a, p, ups1, ...
                                exz, V1, v1, V10, delta, titlestr);
            plot_solution_lambert(vecr1, vecr2, ZKNr, Co, Colors,...
                                  chi1, chi2, a, exz, ups0, titlestr);
        end
    end
end


% -------------------------------------------------------------------------
% Beginn Funktionen
% -------------------------------------------------------------------------                
function ttlprint= PrintOut(ZK, ZKNr, dm, r1, r2,...
                    TOF, gamma, ups0, a, p, ups1,...
                    exz, V1, v1, V10, delta, titlestr)
ttlprint = strjoin(["Lambert-Problem für",titlestr(ZKNr)]);
fprintf('\n %s', ttlprint);
fprintf('\n')
if ZKNr < 3
    fprintf('|  r1 (AE)  |  r2 (AE)  |  TOF (d) | gamma (°) |\n');
    fprintf('|  %7.4f  |  %7.4f  |  %6.1f  | %7.2f   |\n',...
                    r1, r2, TOF, rad2deg(gamma));
    fprintf('\n')
    fprintf('|  a  (AE)  |  p  (AE)  |    e     |  ups0 (°) | ups1 (°) |\n');
    fprintf('|  %7.4f  |  %7.4f  |  %6.4f  |  %7.2f  | %+7.2f  |\n',...
                    a, p, exz, rad2deg(ups0), rad2deg(ups1));
    fprintf('\n')
    DeltaV = (V1-V10);
    v10    = vecnorm(V10);
    deltav = sign(v1-v10)*vecnorm(DeltaV);
    fprintf('| v1 (km/s) | v10 (km/s)| Deltav (km/s)| Winkel(V1,V10)(°)|\n');
    fprintf('|  %7.4f  |  %7.4f  |    %7.4f   |     %7.2f      |\n',...
              v1, v10, deltav,rad2deg(delta));
    fprintf('\n|         V10 (km/s)        |         V1 (km/s)         |      DeltaV (km/s)        |\n');
    fprintf('| [%+7.3f,%+7.3f,%+7.3f] | [%+7.3f,%+7.3f,%+7.3f] | [%+7.3f,%+7.3f,%+7.3f] |\n',...
                    V10, V1, DeltaV);
    fprintf('\n')
else
    fprintf('|  r1 (km)  |  r2 (km)  |  TOF (s) | gamma (°) |\n');
    fprintf('|  %7.1f  |  %7.1f  |  %6.1f  | %7.2f   |\n',...
                    r1, r2, TOF, rad2deg(gamma));
    fprintf('\n')
    fprintf('|  a  (km)  |  p  (km)  |    e     |  ups0 (°) | ups1 (°) |\n');
    fprintf('|  %7.1f  |  %7.1f  |  %6.4f  |  %7.2f  | %+7.2f  |\n',...
                    a, p, exz, rad2deg(ups0), rad2deg(ups1));
    fprintf('\n')
    DeltaV = V1-V10;
    v10    = vecnorm(V10);
    deltav = sign(v1-v10)*vecnorm(DeltaV);
    fprintf('| v1 (km/s) | v10 (km/s)| Deltav (km/s)| Winkel(V1,V10)(°)|\n');
    fprintf('|  %7.4f  |  %7.4f  |    %7.4f   |    %7.2f       |\n',...
                   v1, v10, deltav,rad2deg(delta));
    fprintf('\n|         V10 (km/s)        |         V1 (km/s)         |      DeltaV (km/s)        |\n');
    fprintf('| [%+7.3f,%+7.3f,%+7.3f] | [%+7.3f,%+7.3f,%+7.3f] | [%+7.3f,%+7.3f,%+7.3f] |\n',...
                    V10, V1, DeltaV);
    fprintf('\n')
end
end

% Darstellung Ellipse im Fokuspunkt-KOS
function plot_ell(a, exz, u, uneg, u0, Colors)
  Style = ["-", "-.", ":", ":", "--"];
  r = a*(1-exz^2)./(1+exz*cos(u-u0));
  xe  = r.*cos(u);
  ye  = r.*sin(u);
  plot(xe,ye,'Color', Colors(9,:),'LineWidth',2,'LineStyle',...
      Style(2));
  r = a*(1-exz^2)./(1+exz*cos(uneg-u0));
  xe  = r.*cos(uneg);
  ye  = r.*sin(uneg);
  plot(xe,ye,'Color', Colors(9,:),'LineWidth',1,'LineStyle',...
      Style(3));
end

% Polargleichung Ellipse 
function r = r_ell(p, exz, u, u0)
  r = p./(1+exz*cos(u-u0));
end

function plot_solution_lambert(vecr1, vecr2, ZKNr, Co, Colors,...
                                  chi1, chi2, a, exz, ups0, titlestr)

    figure()
    hold on
    r1 = vecnorm(vecr1); r2 = vecnorm(vecr2); 
    line([0 vecr1(1)], [0 vecr1(2)],...
        'color',Colors(Co(ZKNr,1),:),'LineWidth',2);
    line([0 vecr2(1)], [0 vecr2(2)],...
        'color',Colors(Co(ZKNr,2),:),'LineWidth',2);
    PlotCircleL(0,0,r1,Colors(Co(ZKNr,1),:),1,':');
    PlotCircleL(0,0,r2,Colors(Co(ZKNr,2),:),1,':');
    u    = linspace(chi1,chi2,361);
    uneg = linspace(chi2,chi1+2*pi,361);
    plot_ell(a, exz, u, uneg, ups0, Colors);
    % Apsidenlinie, Periapsis, Apoapsis 
    upsA = pi+ups0;
    xP = a*(1-exz^2)*cos(ups0)/(1+exz*cos(ups0-ups0));
    yP = a*(1-exz^2)*sin(ups0)/(1+exz*cos(ups0-ups0));
    xA = a*(1-exz^2)*cos(upsA)/(1+exz*cos(upsA-ups0));
    yA = a*(1-exz^2)*sin(upsA)/(1+exz*cos(upsA-ups0));
    line([0 xP], [0 yP],'color','k');
    line([0 xA], [0 yA],'color','k');
    % Fokuspunkte
    Foc1=[0 0];
    Foc2=[cos(upsA) sin(upsA)]*2*exz*a;
    PlotCircle(Foc1(1),Foc1(2),0.05*r1,Colors(Co(ZKNr,3),:),2);
    PlotCircle(Foc2(1),Foc2(2),0.02*r1,Colors(15,:),2);
    if ZKNr ==3 
        ylabel('y in km'); xlabel('x in km')
    else
        ylabel('y in AE km'); xlabel('x AE km')      
    end
    % Achsen
    xB = max(r1,r2)*1.1;
    xlim([-xB xB]); ylim([-xB xB]);
    axis equal
    grid on
    ttl=title(strjoin(['Lambert Problem', titlestr(ZKNr)]));
    set(ttl,'FontSize',14,'FontWeight','normal');
    set(gca,'FontSize',16);
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
