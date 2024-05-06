% -------------------------------------------------------------------------
% ReentrySummary.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Dynamik des Wiedereintritts eines Raumschiffes
% in die Erdatmosphäre für verschiedene Eintrittswinkel gammaRe.
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
% Initialisierung
% alle Daten in kg, m und s 
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

gammaRe = deg2rad([2,4]);                     % FPA @ Reentry
vRe     = sqrt((2*muE/rA)*q^2*qminus./(q^2-cos(gammaRe).^2)); % v @ Reentry

% vRe     = vA*vRe./vRe;
Deltav  = cos(gammaRe).*vRe/q - vA;                % Delta v for Reentry
exz     = (q^2-(2*q-1).*(cos(gammaRe)).^2)./(q^2-(cos(gammaRe)).^2);
a       = rA./(1+exz);
rP      = a.*(1-exz);
vARe    = sqrt(muE*(2./rA-1./a));
upsRe   = 180-acosd((qminus^2.*(cos(gammaRe)).^2-q^2.*(sin(gammaRe)).^2)/...
                (qminus^2+q^2.*(sin(gammaRe)).^2));
h0      = hRe;
x0      = 0;
% Parameter der DGL

mL         = 5e03;                % Payload in kg
rho0       = 1.752e09;            % Dichte in kg/km^3
H0         = 6.7;                 % H0 in km
kappaD     = 25;                  % normierte Widerstandswert
cwD        = 2*kappaD*mL/rho0/H0; % Drag Widerstandsbeiwert x Fläche km^2 
facAD      = 0.3;                 % Verhältnis Auftrieb/Drag
cwAvec     = facAD*cwD;           % Auftriebswiderstandswert x Fläche km^2 

% Barometrische Höhenformel
rho        = @(h)rho0*exp(-h/H0);  

% AB für DGL
P1.mL     = mL;
P1.gE     = gE;
P1.RE     = RE;
P1.cwD    = cwD;
P1.facAD  = facAD;

% MATLABs Runge-Kutta ode45 Routine 
opts   = odeset('AbsTol',1.e-9,'RelTol',1.e-7);
% Schleife über verschieden Eintrittswinkel gammaRe
for k=1:length(gammaRe)
    tmax(k)  = round(1.0*sqrt(2*hRe/gE/sin(gammaRe(k)))/100)*100;
    xmax(k)  = round(1.0*hRe/tan(gammaRe(k))/200)*200;
    hmax = 125;
    tspan1 = linspace(0,tmax(k),1001);
    AB1    = [vRe(k);gammaRe(k);h0;x0];          % AB für DGL
    [t,Y1]      = ode45(@DGL_Reentry, tspan1, AB1, opts, P1, rho);
    v1(k,:)     = Y1(:,1);
    gamma1(k,:) = Y1(:,2);
    x1(k,:)     = Y1(:,4);
    h1(k,:)     = Y1(:,3);
    t1(k,:)     = t(:);
end  
% Berechnung Beschleunigungen
for k=1:length(gammaRe)
    a1(k,:) = zeros(1,length(v1(k,:)));
    for m=2:length(v1(k,:))
        a1(k,m) = (v1(k,m)-v1(k,m-1))/(t1(k,m)-t1(k,m-1))/gE;
    end
    a1(k,1) = a1(k,2);
end

%% Graphische Ausgabe

strkappa = string(num2str(kappaD,'%-2.0f'));
strcDA   = string(num2str(cwD,'%4.1e'));
strfacAD = string(num2str(facAD,'%3.1f'));
strPara  = strjoin(['Parameter: kD=',strkappa,...
              ';  cwDA:=',strcDA, 'km^2;  A/D=',strfacAD']);

%----------------------------------------------------------------------
figure('name','h-x-Plot')
hold on
for k=1:length(gammaRe)
    p(k)=plot(x1(k,:),h1(k,:),'color',Colors(3,:),...
         'Linewidth',2,'LineStyle',Style(2*(k-1)+1));
    lgdstr(k,:) = strjoin(["h(x) für \gamma_{Re} = ",...
                  string(num2str(rad2deg(gammaRe(k)),'%3.0f')),'°']);
end    
grid on
grid on;
text(max(xmax)/5,hRe-10,strPara,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p,lgdstr,'Location','northeast','NumColumns',1); 
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, max(xmax), 0 hmax]);
xlabel('x in km','FontSize',14); 
ylabel('h in km','FontSize',14)
legend box off;
set(gca,'FontSize',14);

%----------------------------------------------------------------------
figure('name','Zeit-Plot')
%++++++++++++++++++++++++++
subplot(2,1,1)
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, max(tmax)-100, 0 hRe]);
hold on
for k=1:length(gammaRe)
    pl(k)= plot(t1(k,:),h1(k,:),'color',Colors(3,:),...
           'Linewidth',1,'LineStyle',Style(2*(k-1)+1));
    lgdstrleft(k,:) = strjoin(["h(t) für \gamma_{Re} = ",...
                  string(num2str(rad2deg(gammaRe(k)),'%3.0f')),'°']);
end 
text(max(tmax)/5,hRe-10,strPara,'FontSize',14,'FontWeight','normal'); 
hp1=legend(pl,lgdstrleft,'Location','northeast','NumColumns',1); 
legend box off;
set(hp2,'FontSize',14,'FontWeight','normal'); 
xlabel('t in s','FontSize',14); 
ylabel('h in km','FontSize',14)
grid on
grid minor
set(gca,'FontSize',16);
%++++++++++++++++++++++++++
subplot(2,1,2);
% .......................................................
yyaxis left
set(hp2,'FontSize',14,'FontWeight','normal'); 
axis([0, max(tmax)-100, 0 1.2*max(vRe)]);
hold on
for k=1:length(gammaRe)
    p2(k)= plot(t1(k,:),v1(k,:),'color',Colors(3,:),...
           'Linewidth',1,'LineStyle',Style(2*(k-1)+1));
   lgdstr(k,:) = strjoin(["v(t) für \gamma_{Re} = ",...
                  string(num2str(rad2deg(gammaRe(k)),'%3.0f')),'°']);
end 
ylabel('v in km/s','FontSize',14)
set(gca,'FontSize',14);
% .......................................................
yyaxis right
hold on
axis([0 max(tmax)-100, 0 round(1.2*max(max(abs(a1))))]);
xlabel('t in s','FontSize',14); 
ylabel('|a| in g','FontSize',14)
for k=1:length(gammaRe)
    p2(k+length(gammaRe)) = plot(t1(k,:),-a1(k,:),'color',Colors(5,:),...
           'Linewidth',1,'LineStyle',Style(2*(k-1)+1));
    lgdstr(k+length(gammaRe),:) = strjoin(["a(t) für \gamma_{Re} = ",...
                  string(num2str(rad2deg(gammaRe(k)),'%3.0f')),'°']);
end    
hp3=legend(p2, lgdstr,'Location','northeast','NumColumns',1); 
set(hp3,'FontSize',14,'FontWeight','normal'); 
legend box off;
grid on
grid minor
set(gca,'FontSize',14);


%%

% DGL
function dY = DGL_Reentry(t, Y, P1, Dichte)
% Y(1):v, Y(2):gamma, Y(3):h, Y(4):x 
mR    = P1.mL;
F_D   = 0.5*feval(Dichte,Y(3))*P1.cwD*Y(1)^2;
F_A   = F_D*P1.facAD;
g     = P1.gE./(1+2*Y(3)./P1.RE+(Y(3)/P1.RE)^2);
dY    = [-F_D/mR + g*sin(Y(2));
         (-F_A/mR/Y(1) +(g/Y(1)-Y(1)/(P1.RE+Y(3)))*cos(Y(2)));
         - Y(1)*sin(Y(2));
         Y(1)*cos(Y(2))];
end

