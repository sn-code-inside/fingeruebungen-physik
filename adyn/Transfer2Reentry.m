% -------------------------------------------------------------------------
% Transfer2Reentry.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet den Transferbahnparameter aus einer
% Kreisbahn zum Wiedereintritts eines Raumschiffes in 
% die Erdatmosphäre in Abhängigkeit verschiedener Eintrittswinkel gammaRe.
%
%--------------------------------------------------------------------------


%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%%
% Initialisierung % alle Daten in kg, m und s 
G   = 6.67430e-11;                 % G in m3 /kg /s2
ME  = 5.98e24;                     % Masse der Erde in kg
RE  = 6.378e3;                     % Erdradiuse in km
muE = G*ME*1e-9;                   % µ in km3 /kg /s2
gE  = muE/RE^2;                    % gE in km/s^2

hA  = 422;                         % Anfangshöhe in km
hRe = 122;                         % Reentry Höhe in km
rA      = RE + hA;
rRe     = RE + hRe;
vA      = sqrt(muE/rA);
q       = rA/rRe;
qminus  = q-1;

%% Anfangswerte 
prompt   = '\nBitte Flight Path Angle in ° eingeben oder [E] für Exit :  ';
gammaRestr = input(prompt,'s');
gammaRe  = str2double(gammaRestr);
x = ~isnan(gammaRe); 
if x==1
    gammaRe = str2double(gammaRestr);
    gammaRe = deg2rad(gammaRe); 
    cosgam  = cos(gammaRe);
    singam  = sin(gammaRe);
    % Exzentrizität, a rP
    exz     = (q^2-(2*q-1).*cosgam.^2)./(q^2-cosgam.^2);
    a       = rA./(1+exz);
    p       = a*(1-exz^2);
    rP      = a.*(1-exz);
    % v @ Reentry
    vRe     = sqrt((2*muE/rA)*q^2*qminus./(q^2-cos(gammaRe).^2));
    % Delta v for Reentry
    Deltav  = cosgam.*vRe/q - vA;                
    vARe    = sqrt(muE*(2./rA-1./a));
    cosupsW = (qminus^2*cosgam.^2-q^2*singam.^2)./(qminus^2+(2*q-1)*singam.^2);
    cosups  = (q*(1-exz)-1)./exz;
    thetRe  = abs(acos(cosups)-pi);     
    
    %% Graphische Ausgabe

    strhA    = string(num2str(hA,'%4.0f km'));
    strPara  = strjoin(['Parmeter: h_A=',strhA,...
                    ' \gamma_{Re}= ',gammaRestr, '°']);
    u     = linspace(0,2*pi,3601);
    sinu  = sin(u);
    cosu  = cos(u);
    xA    = rA*cosu/RE;
    xE    = cosu;
    yA    = rA*sinu/RE;
    xRe   = rRe*cosu/RE;
    xE    = cosu;
    yRe   = rRe*sinu/RE;
    yE    = sinu;
    xTr   = p*cosu./(1+exz*cos(u-1*pi))/RE;
    yTr   = p*sinu./(1+exz*cos(u+1*pi))/RE;
    xTrRe = p*cos(thetRe)./(1+exz*cos(thetRe-1*pi))/RE;
    yTrRe = p*sin(thetRe)./(1+exz*cos(thetRe-1*pi))/RE;
    maxAx = 1.1*rA/RE;          
    %----------------------------------------------------------------------
    figure(1)
    hold on
    hp(1)=plot(xA,yA,'color',Colors(2,:),'Linewidth',1,'LineStyle',Style(1));
    hp(2)=plot(xE,yE,'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(1));
    hp(3)=plot(xTr,yTr,'color',Colors(4,:),'Linewidth',2,'LineStyle',Style(1));
    hp(4)=plot(xRe,yRe,'color',Colors(8,:),'Linewidth',1,'LineStyle',Style(3));
    line([0 rA/RE],[0 0]);  
    line([0 xTrRe],[0 yTrRe]);
    PlotCircle (0,0,0.01,Colors(3,:),1);
    axis equal
    axis(maxAx*[-1.1 1.1, -1.1 1.1]);
    xlabel('x in RE','FontSize',14); 
    ylabel('y in RE','FontSize',14);
    grid on
    title(strPara,'FontSize',12,'FontWeight','normal');
    hp2=legend(hp,'Anfangsorbit (ISS)','Erde', 'Transferbahn (Deorbit)',...
               'Reentry-Höhe','Location','bestoutside','NumColumns',1); 
    set(hp2,'FontSize',12,'FontWeight','normal'); 
    legend box off;
    set(gca,'FontSize',12);

% PrintOut
ttlprint = "Transfer zum Wiedereintritt (Deorbit)";
fprintf('\n %s \n', ttlprint);
fprintf('\n |   hA  (km) |   rA  (km) |   q=rA/rRE   | gamma_Re (°) |   vA  (km/s) |  \n');
fprintf(  ' |  %7.1f   |  %7.1f   |   %7.3f    |   %7.3f    |   %7.3f    |',...
       hA, rA, q,  rad2deg(gammaRe), vA);
fprintf('\n | vRe (km/s) | theta_Re (°) | Deltav (km/s) |  Exz. e  |\n');
fprintf(  ' |  %7.3f   |   %7.2f    |   %7.3f     |  %6.4f  | \n',...
       vRe, rad2deg(thetRe), Deltav, exz);
fprintf('\n');

end




% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------