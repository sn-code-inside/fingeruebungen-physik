% -------------------------------------------------------------------------
% TransferUnendlich.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Geschwindigkeitsbedarfe 
% für Hohmannn-Oberth-Edelbaum Transfers nach Unendlich, ebenso wie die 
% FLugzeiten (normiert auf die Anfangsgeschwindigkeit/Umlaufzeit/Radius
% der Anfangs-Kreisbahn)
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


% Parameter
GME = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius in km


%% Berechnung nach Geschwindigkeitsbedarfs

vinfr = linspace(0,3,91);  %Verhältnis vinf/v1
 

% Hohmann-Transfer
delv_H = sqrt(2 + vinfr.^2)-1;

% Oberth-Transfer
r1min = 20;  %  = r0/r_
delv_O = 1+ sqrt(vinfr.^2 + 2*r1min)-sqrt(2+2*r1min);

% Edelbaum-Transfer
r1max = 0.4; %  = r0/r+
delv_E = (sqrt(vinfr.^2 + 2*r1min)+ sqrt(2+2*r1max)-sqrt(2*r1min+2*r1max)-1);

% Graphik Geschwindigkeitsbedarf

figure(1)
plot(vinfr,delv_H, 'linewidth',2,'color',Colors(2,:));
hold on
plot(vinfr,delv_O, 'linewidth',2,'color',Colors(3,:));
plot(vinfr,delv_E, 'linewidth',2,'color',Colors(4,:));
grid on
ylabel('\Delta v / v_1')
xlabel('v_{inf} / v_1')
legend('Hohmann','Oberth','Edelbaum',...
        'location','northwest');   
legend box off
ttl=title('Transfers nach Unendlich');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',16);

%% Berechnung der Flugzeiten

% Kreisbahn initial
r1 = 7000;
v1 = sqrt(GME/r1);
T1 = 2*pi*r1/v1;
T1min = T1/60;

% Variation finaler Abstand
rEnd1 = 50*r1;  
rEnd2 = 100*r1;
rEnd3 = 200*r1;
delv = 1.1*v1;  %Geschwindigkeitsbudget

% Ellipsenparameter
rOut = r1*linspace(1,10,200);   %Apogäum der Edelbaum-Ellipse
rIn  = 0.05*r1;                 %Perigäum der Oberth-Ellipse


% Hohmann-Transfer
rin   = r1;
rout  = r1*rOut./rOut;
vEnd  = vfinal(v1, delv, rin, rout, r1); % Berechnung Endgegeschwindigkeit 
vinfR = vEnd/v1; 
time_H1 = thyp(vinfR, rEnd1, r1, r1);
time_H3 = thyp(vinfR, rEnd3, r1, r1);

% Oberth-Transfer
rin   = rIn;
rout  = r1*rOut./rOut;
vEnd  = vfinal(v1, delv, rin, rout, r1); % Berechnung Endgegeschwindigkeit 
vinfR = vEnd/v1; 
time_O1 = thyp(vinfR, rEnd1, r1, rin) + tell(r1,r1,rin);
time_O3 = thyp(vinfR, rEnd3, r1, rin) + tell(r1,r1,rin);

% Edelbaum-Transfer
rin   = rIn;
rout  = rOut;
vEnd  = vfinal(v1, delv, rin, rout,r1); % Berechnung Endgegeschwindigkeit 
vinfR = vEnd/v1; 
time_E1 = thyp(vinfR, rEnd1, r1, rin) + tell(r1,r1,rout)+ tell(r1,rout,rin);
time_E3 = thyp(vinfR, rEnd3, r1, rin) + tell(r1,r1,rout)+ tell(r1,rout,rin);

% Graphik Flugzeiten

figure(2)
subplot(1,2,1)
plot(rOut/r1,time_H1,'linewidth',2,'color',Colors(2,:));
hold on
plot(rOut/r1,time_O1,'linewidth',2,'color',Colors(3,:));
plot(rOut/r1,time_E1,'linewidth',2,'color',Colors(4,:));
grid on
ylabel('t_F / T_1')
xlabel('r_+ / r_1')
axis([0 max(rout/r1) min(time_E1)*0.9 max(max(time_H3),max(time_E3))*1.1]); 
legend('Hohmann','Oberth','Edelbaum',...
        'location','northwest');   
legend box off
ttl=title('r_{end} = 50 r_1');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);

subplot(1,2,2)
plot(rOut/r1,time_H3,'linewidth',2,'color',Colors(2,:));
hold on
plot(rOut/r1,time_O3,'linewidth',2,'color',Colors(3,:));
plot(rOut/r1,time_E3,'linewidth',2,'color',Colors(4,:));
grid on
ylabel('t_F / T_1')
xlabel('r_+ / r_1')
axis([0 max(rout/r1) min(time_E1)*0.9 max(max(time_H3),max(time_E3))*1.1]); 
legend('Hohmann','Oberth','Edelbaum',...
        'location','southwest');   
legend box off
ttl=title('r_{end} = 200 r_1');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);


% ------------------------------------------------------------------------% 
%-------------------------------------------------------------------------%
% Ende Programm
% ------------------------------------------------------------------------%


%% Funktionen
% Flugzeit auf der Hyperbel
function th = thyp(vinf, rinf, r1, rP)
    fac1 = 1/2/pi./vinf.^3;
    sum1 = sqrt((1+vinf.^2.*rinf/r1).^2 - (1+vinf.^2.*rP/r1).^2);
    sum2 = acosh((1+vinf.*rinf/r1)./(1+vinf*rP/r1));
    th = fac1.*(sum1-sum2);
end
% Flugzeit auf der halben Ellipse
function te = tell(r1,rA,rP)
    te = sqrt(((rA+rP)/2/r1).^3)/2;
end

% Endgeschwindigkeit
function ve = vfinal(v1, delv, rin, rout, r1)
    ve = v1*sqrt( (1+delv/v1 + sqrt(2*r1/rin+2*r1./rout)-...
                       sqrt(2+2*r1./rout)).^2 - 2*r1/rin);
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
