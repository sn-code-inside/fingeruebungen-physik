% -------------------------------------------------------------------------
% Raketenstart01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet Start einer einstufigen Rakete von der Erdoberfläche
%
%--------------------------------------------------------------------------


%%
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;

%%
% Initialisierung
% alle Daten in kg, m und s 
ME = 5.98e24;                     % Masse der Erde in kg
RE = 6.378e6;                     % Erdradiuse in m
G  = 6.671e-11;                   % G in (m^3 / kg /s^2)
P1.mE = ME;                       % Erdmasse in kg
P1.G  = G;

% Normierung
tN    = sqrt(RE^3/G/ME);          % Umlaufzeit Erdhöhe, Einheit Zeit in s
vN    = RE/tN;                    % Einheit der Geschwindigkeit               
FN    = ME*RE/tN^2;               % Einheit der Kraft               

mi = 2.8e6;                       % Anfänglicher Payload+Treibstoff in kg 
fT = 0.96;  mT = mi*fT;           % Anteil und Masse Treibstoff (norm.)

Schub = 1.5*(G*mi*ME/RE^2);       % Schub = 1.5 x initiales Gewicht in N
u     = 3000;                     % Austrittsgeschwindigkeit m/s 

theta  = deg2rad(55);             % Start Winkel = Brennwinkel 
P1.ux=u*cos(theta); 
P1.uy=u*sin(theta);               

alpha = Schub/u;                  % Treibstoffverbrauch 
P1.alpha = alpha;

mf     = mi-mT;                   % Finale Masse = Payload (nach burnout)
tf     = (mi-mf)/alpha;           % Treibstoff burnout Zeit in
P1.tf  = tf;

tend   = 500*tf;
tspan  = linspace(0,tend);        % Simulationsdauer

x0  = 0; y0  = RE; vx0 = 0; vy0 = 0;                
AB = [x0;vx0;y0;vy0;mi];          % AB für DGL

% MATLABs Runge-Kutta ode45 Routine 
opts = odeset('AbsTol',1.e-7,'RelTol',1.e-5);
[t,Y]=ode45(@(t,Y, P1)DGL_Rakete(t,Y,P1),[0 tend],AB,opts,P1);


%%
% Graphische Ausgabe
xmax = max(Y(:,1));
ymax = max(Y(:,3));

n=length(t);                  
for i=1:n                          % Auswahl Punkte ober Erdoberfläche
  if sqrt(Y(i,1)^2+Y(i,3)^2) >= 0.99*RE 
     nn=i;
     t1(i)=t(i);
     x1(i) = Y(i,1); y1(i) = Y(i,3);
     vx1(i)= Y(i,2); vy1(i)= Y(i,4);
  else
      break;
  end
end
tplotmax = t(nn);

% Plot Geschwindigkeit über Zeit
v1=sqrt(vx1.^2+vy1.^2);
r1=sqrt(x1.^2+y1.^2);
subplot(2,1,1)
plot(t1,v1/1000,'Color',Colors(2,:),'Linewidth',2); hold on
plot(t1,vx1/1000,'Color',Colors(3,:),'Linewidth',2);
plot(t1,vy1/1000,'Color',Colors(4,:),'Linewidth',2);
ylabel('v, v_x, v_y in km/s','FontSize',14)
str= "Geschwindigkeit in km/s";
grid on;
h1 = title(str,'FontSize',12);
set(h1,'FontSize',14,'FontWeight','normal'); 
h2=legend('v','v_x','v_y'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
axis([0, 1.2*tplotmax, 1.2*min(min(vy1),min(vx1))/1000 1.2*max(v1)/1000]);
grid on;
legend box off;
set(gca,'FontSize',16);

% Plot Abstände über Zeit
subplot(2,1,2)
% plot(t1,(r1-RE)/1000,'Color',Colors(2,:),'Linewidth',2); hold on
plot(t1,r1/1000,'Color',Colors(2,:),'Linewidth',2); hold on
plot(t1,x1/1000,'Color',Colors(3,:),'Linewidth',2);
plot(t1,y1/1000,'Color',Colors(4,:),'Linewidth',2);
xlabel('t in s','FontSize',14); ylabel('r, x, y in km','FontSize',14)
h2=legend('r','x','y'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
str= 'Entfernung vom Erdmittelpunkt in km';
axis([0, 1.2*tplotmax, 1.2*min(min(x1),min(y1))/1000, 1.2*max(r1)/1000]);

% str= 'Entfernung von der Erdoberfläche in km';
% axis([0, 1.2*tplotmax, 0, 1.2*(max(r1)-RE)/1000]);
h1 = title(str,'FontSize',12);
set(h1,'FontSize',14,'FontWeight','normal'); 

grid on;
legend box off;
set(gca,'FontSize',16);

%% 
% Flugsimulation 
figure()
hold on
L=1.25*sqrt(xmax^2+ymax^2)/1000;      %Plotgröße
phi = [0:0.025:2*pi];
xE=RE*cos(phi);yE=RE*sin(phi);             %Erdumfang
axis ([-L L -L L])           %windows size
axis square                  %square window
xlabel('x in km ','FontSize',14);
ylabel('y in km ','FontSize',14);
str= ' Simulation ';
h1 = title(str,'FontSize',12);
set(h1,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);
grid on;
str2=cat(2,' Burnout-Zeit = ',num2str(tf,3),...
           ' s ,  Startwinkel = ',num2str(rad2deg(theta),3),'^o');
str3=cat(2,' Schub = ',num2str(Schub/1000,5),' kN,  Treibgasgeschw. = ',...
            num2str(u/1000,3),' km/s');
text(-0.95*L,+L*(1-0.125),str2,'FontSize',12)
text(-0.95*L,+L*(1-0.25),str3,'FontSize',12)
warning off;
Range=RE*abs(atan(y1(nn)/x1(nn))-atan(y0/x0));
str6=cat(2,' \alpha = ',num2str(alpha,3),...
           ' kg/s , Reichweite = ',num2str(Range,3),' m');
text(-0.95*L,+L*(1-0.375),str6,'FontSize',12)
hp(1) = plot(xE/1000,yE/1000,'color',Colors(3,:),'Linewidth',2); %Erde
hp(2) = plot(x1/1000,y1/1000,'color', Colors(4,:),...
        'Linewidth',2, 'LineStyle', ':'); %Trajektorie
i=1;
hp(3) = plot(x1(i)/1000,y1(i)/1000,'d','color', Colors(4,:)); % Rakete
h = hp(3);
% for i=2:nn                     
%   h.Visible = 'off';
%   h = plot(x1(i)/1000,y1(i)/1000,'d','color', Colors(4,:));%Position Rakete
%   h.Visible = 'on';
%   pause(0.1)
% end
% h = legend(hp,'Erde','Trajektorie', 'Rakete','location','south'); 
% set(h,'FontSize',12); legend box off;


%%

% DGL
function dY = DGL_Rakete( t, Y, P1)
% Y(1):x, Y(2):vx, Y(3):y, Y(4):vy, Y(5):m
r12=sqrt(Y(1).^2+Y(3).^2);
tmp=P1.alpha*stepf(t,P1.tf,1);
dY = [Y(2);P1.ux*tmp./Y(5)-P1.G*P1.mE*Y(1)./r12.^3; Y(4);...
               P1.uy*tmp./Y(5)-P1.G*P1.mE*Y(3)./r12.^3;-tmp;];
end

