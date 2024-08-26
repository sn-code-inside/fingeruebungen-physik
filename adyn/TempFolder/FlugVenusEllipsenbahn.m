% -------------------------------------------------------------------------
% FlugVenusEllipsenbahn.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Ellipsen-Flugbahn einer Raumsonde zu einem der 
% inneren Planeten (Venus, Merkur) zu verschiedenen Zeiten
% auf Basis der der Lösung des Lambert-Problems.
% Vorgegeben sind das Startdatum und das geplante Ankunftsdatum.
% Gesucht ist die Lambert-Transfer-Bahn und das notwendige Delta v.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Initialisierung

% Bitte ein Datum auswählen
% dt1 = datetime('2022-06-18 00:00:00');
% dt2 = datetime('2022-08-10 00:00:00');
dt1 = datetime('2023-09-14 00:00:00');
dt2 = datetime('2024-02-20 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum %Startdatum
DatumStr1 =  string(dt1,'dd.MM.yy');
T2 = juliandate(dt2); % Julianisches Datum %Startdatum
DatumStr2 =  string(dt2,'dd.MM.yy');
TOF       = T2-T1; % Flugzeit in Tagen
Aequi = 'Datum';
% Einlesen der Bahnparameter
[BaPa,BaPadot]=OrbitParameter(T1,Aequi);
GS  = 2.95479E-04;     % µG [AE^3/Tage^2] für Sonne
G   = 6.671e-11;       % G in in m^3/s^2/kg
MS  = 1.989e30;        % kg
muSk= G*MS*1e-9;       % in km^3/s^2

%% Vorbereitung Rechnung

maxT    = 365; % Berechnungsbereich
scFac   = 4;   % Schrittweite 6 h 
T_Fv    = linspace(T1,T2,maxT*scFac+1); %Flugzeitvektor 
maxT    = 365; % Berechnungsbereich Jahr
scFac   = 4;   % Schrittweite 6 h 
T_Yv  = linspace(T1,T1+maxT,maxT*scFac+1); %Jahr Erde
calendardays = linspace(0,maxT,maxT*scFac+1);
ianf = 2;
iend = 3;

% Berechnung nach Keplerloesung
for k=1:4
    PlanetsY(k)=PlanetPQR(T_Yv, BaPa, BaPadot, k);
    PlanetsF(k)=PlanetPQR(T_Fv, BaPa, BaPadot, k);
end

% Bestimmung der Parameter für das Lambert-Problem
% Wir vereinfachen uns das Problem, in  dem wir die Bewegung in der 
% Ekliptik-Ebene annehmen
rE1   = PlanetsF(3).ekl(1,1);
upsE1 = rad2deg(PlanetsF(3).ekl(2,1));
rV2   = PlanetsF(2).ekl(1,end);
upsV2 = wrapTo360(rad2deg(PlanetsF(2).ekl(2,end)));
rV1   = PlanetsF(2).ekl(1,1);
upsV1 = rad2deg(PlanetsF(2).ekl(2,1));
rE2   = PlanetsF(3).ekl(1,end);
upsE2 = wrapTo360(rad2deg(PlanetsF(3).ekl(2,end)));


gammad  = winkel(PlanetsF(2).xyz(:,end),PlanetsF(3).xyz(:,1));


%% Rechnung Lambert, Iteration -------------------------------------------% 

r1    = rE1;          % in AE
r2    = rV2;          % in AE
gamma = upsV2-upsE1;

%Iteration, Bisektion
[aiter, amin, amax, alpha_n, beta_n, tau_n, tm, tpar, nend] = ...
    BiSektLamb1(r1, r2, gamma, TOF, GS);

% PrintOut
if TOF-1 < tpar 
   fprintf('\n  Kein elliptischer Orbit möglich, \n');
   fprintf('\n  da t_min > gegebener TOF. \n');          
else
% Berechnung Transferbahn (Halbparameter p, Elliptizität, Perihelrichtung)
[xtr, ytr, p_n, exzf, af, pf, ups0]  = lambert_bahn(aiter, nend, alpha_n,...
                                  beta_n, rE1, rV2, gamma,...
                                  upsE1, upsV2, GS, muSk);
fprintf('\n    n    |  a_n  (AE)  |  p_n  (AE)  |   tau (d)    |   e_n     |\n');
for k=1:nend
    psi = (alpha_n(k)-beta_n(k))/2;
    phi = (alpha_n(k)+beta_n(k))/2;
    fprintf('  %3u    |  %7.5f    |  %7.5f    | %9.2f    |  %6.4f   |\n',...
                k, aiter(k), p_n(k), tau_n(k+1),exzf);
end
end




%% Graphische Ausgabe
%-------------------------------------------------------------------------%

% Trajektorien
header1='Innere Planeten: Projektion Bahnen auf die Ekliptik';
figure('Name',header1);
%Planetenbahnen
for iPlot = ianf:iend
    plot(PlanetsY(iPlot).xyz(1,:),PlanetsY(iPlot).xyz(2,:),'Color', ...
         Colors(iPlot,:));
    hold on
    axis equal
end
% Sonne 
PlotCircle(0,0,0.07,Colors(10,:),2);  
%Positionen Start und Ankunft
k = 1;
for iPlot =  ianf:iend
    p(k)=plot(PlanetsY(iPlot).xyz(1,1),PlanetsY(iPlot).xyz(2,1),'+-',...
              'Color', Colors(iPlot,:),'LineWidth',2,'MarkerSize',8,...
              'MarkerFaceColor',Colors(iPlot,:));
    k=k+1;
    plot([0 PlanetsF(iPlot).xyz(1,1)],[0 PlanetsF(iPlot).xyz(2,1)],...
         ':','Color', Colors(iPlot,:),'LineWidth',1);
    p(k)=plot(PlanetsF(iPlot).xyz(1,end),PlanetsF(iPlot).xyz(2,end),'d-',...
              'Color', Colors(iPlot,:),'LineWidth',2,'MarkerSize',8,...
              'MarkerFaceColor',Colors(iPlot,:));
    k=k+1;
    plot([0 PlanetsF(iPlot).xyz(1,end)],...
         [0 PlanetsF(iPlot).xyz(2,end)],...
          ':','Color', Colors(iPlot,:),'LineWidth',1);
end
k=k+1;
p(5)=plot(xtr,ytr,'Color',...
           Colors(1,:),'LineWidth',2); 
ylim([-1.1 1.1])
xlim([-1.1 1.1]);
grid on;
grid minor,
header2 = strjoin([DatumStr1," - ", DatumStr2]);
ttl=title(header2);
set(ttl,'FontSize',14, 'FontWeight','normal');
xlabel('x in AE')
ylabel('y in AE');
legend(p, strjoin([PlanetsY(ianf).Name, "(bei Start)"]),...
          strjoin([PlanetsY(ianf).Name, "(bei Ankunft)"]),...
          strjoin([PlanetsY(iend).Name, "(bei Start)"]),...
          strjoin([PlanetsY(iend).Name, "(bei Ankunft)"]),'Transferbahn',...
          'location','bestoutside','numcolumns',1);
                %"Transferbahn",

legend boxoff;
set(gca,'FontSize',16);


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

function ups = winkel(vec1,vec2)
    ups = acosd(dot(vec1,vec2)/vecnorm(vec1)/vecnorm(vec2));
end
 
% Iteration/Bisektion
function [aiter,amin,amax,alpha_n,beta_n,tau_n,tm,tpar,nend] = ...
          BiSektLamb1(r1, r2, gamma, TOF, GZK)
    c     = sqrt(r1^2+r2^2-2*r1*r2*cosd(gamma));
    s     = (c+r1+r2)/2;
    if gamma < 180
        fac_beta  = 1;
    else
        fac_beta = -1;
    end
    amin0    = s/2;
    alpham  = wrapTo2Pi(2*asin(sqrt(s/2/amin0)));
    alphamd = rad2deg(alpham);
    betam   = fac_beta*2*asin(sqrt((s-c)/2/amin0));
    betamd  = rad2deg(betam);
    tm      = sqrt(amin0^3/GZK)*(alpham-betam-(sin(alpham)-sin(betam)));
    tpar    = (sqrt(2)/3)*sqrt(s^3/GZK)*(1-fac_beta*((s-c)/s)^(3/2));
    if TOF < tm
        fac_alpha = 1;
        sum_alpha = 0;
    else
        fac_alpha = -1;
        sum_alpha = 2*pi;
    end
    amax = zeros(1,100);
    amin = zeros(1,100);
    amax(1)    =  5*amin0;
    amin(1)    =  1*amin0;
    aiter   = zeros(1,100);
    beta_n  = zeros(1,100);
    alpha_n = zeros(1,100);
    tau_n   = zeros(1,100);
    aiter(1)= (amin(1)+amax(1))/2;
    n       = 1;
    nend    = 100;
    while abs(1-tau_n(n)/TOF) > 0.0001 && n <100
        alpha_n(n)  = wrapTo2Pi(sum_alpha+fac_alpha*2*asin(sqrt(s/2/aiter(n))));
        beta_n(n)   = fac_beta*2*asin(sqrt((s-c)/2/aiter(n)));
        tau_n(n+1)  = sqrt(aiter(n)^3/GZK)*...
                   (alpha_n(n)-beta_n(n)-(sin(alpha_n(n))-sin(beta_n(n))));
        if tau_n(n+1) - TOF  < 0 
            amax(n+1) = amax(n);
            amin(n+1) = aiter(n);
        else
            amax(n+1) = aiter(n);
            amin(n+1) = amin(n);
        end
        aiter(n+1)= (amin(n+1)+amax(n+1))/2;
        n = n+1;
    end
    nend  = n-1;
    if nend >98
        amax = zeros(1,100);
        amin = zeros(1,100);
        amax(1)    =  1*amin0;
        amin(1)    =  5*amin0;
        aiter   = zeros(1,100);
        beta_n  = zeros(1,100);
        alpha_n = zeros(1,100);
        tau_n   = zeros(1,100);
        aiter(1)= (amin(1)+amax(1))/2;
        n       = 1;
        nend    = 100;
        while abs(1-tau_n(n)/TOF) > 0.0001 && n <100
            alpha_n(n)  = wrapTo2Pi(sum_alpha+fac_alpha*2*asin(sqrt(s/2/aiter(n))));
            beta_n(n)   = fac_beta*2*asin(sqrt((s-c)/2/aiter(n)));
            tau_n(n+1)  = sqrt(aiter(n)^3/GZK)*...
                       (alpha_n(n)-beta_n(n)-(sin(alpha_n(n))-sin(beta_n(n))));
            if tau_n(n+1) - TOF  < 0 
                amax(n+1) = amax(n);
                amin(n+1) = aiter(n);
            else
                amax(n+1) = aiter(n);
                amin(n+1) = amin(n);
            end
            aiter(n+1)= (amin(n+1)+amax(n+1))/2;
            n = n+1;
        end
        nend  = n;
    end
end
    
function [xtr, ytr, p_n, exzf, af, pf, ups0]  = lambert_bahn(aiter, nend, ...
                                           alpha_n, beta_n, rE1, rV2, ...
                                           gamma, upsE1, upsV2, GS, muSk)
    af     = aiter(nend);
    alphaf = alpha_n(nend);
    betaf  = beta_n(nend);
    c      = sqrt(rE1^2+rV2^2-2*rE1*rV2*cosd(gamma));
    s      = (c+rE1+rV2)/2;
    p_n    = zeros(1,100);
    if gamma <180
       p_n = (2*aiter*rE1*rV2/c^2)*(1-cosd(gamma)).*(sin((alpha_n+beta_n)/2)).^2;
    else
       p_n = (2*aiter*rE1*rV2/c^2)*(1-cosd(gamma)).*(sin((alpha_n-beta_n)/2)).^2;  
    end
    pf = p_n(nend);
    vec_r1   = [rE1*cosd(upsE1);rE1*sind(upsE1);0]; % in km 
    vec_r2   = [rV2*cosd(upsV2);rV2*sind(upsV2);0]; % in km 
    vec_vtr1 = sqrt(GS*pf)/rE1/rV2/sind(gamma)*...
               ((vec_r2-vec_r1)+rV2/pf*(1-cosd(gamma)*vec_r1));
    % Achtung alles in km und s
    vec_e1 = (1/rE1-1/af)*vec_r1 -(1/muSk)*dot(vec_r1,vec_vtr1)*vec_vtr1;
    exzf   = vecnorm(vec_e1);
    E1    = acosd((af-rE1)/af/exzf);
    ups0  = upsE1 - 2*atand(sqrt((1+exzf)/(1-exzf)*tand(E1/2)));
    usp0  = wrapTo360(ups0);
    upstr = linspace(upsE1,upsE1+360,1000);
    rtr   = af*(1-exzf^2)./(1+exzf*cosd(upstr-ups0));
    xtr   = rtr.*cosd(upstr);
    ytr   = rtr.*sind(upstr);
end 


% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
