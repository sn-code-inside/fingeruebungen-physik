% -------------------------------------------------------------------------
% OberthSwingBy.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Swing-By Manöver mit Oberth Effekt
% am Beispiel Jupiter.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
LW = 'linewidth';
LS = 'linestyle';
LC = 'color';

%% Parameter
% Variablen, Konstanten, Parameter, hier alles in m, kg, s
% Parameter hier alles in m, kg, s
G   = 6.671e-11;      % G in in m^3/s^2/kg
AE  = 149597870700;   % in m
AEk = AE/1000;        % in km
mJ  = 1898.7*1e24;    % Jupitermasse in kg
MS  = 1.989e30;       % Sonnenmasse in kg

% Planetendaten ab hier alles in km,s,kg
muS  = MS*G*1e-9;               % mu Sonne in km^3/s^2
muJ  = 1898.7*1e24*G*1e-9;      % mu Jupiter in km^3/s^2
aJ   = 5.20336301*AEk;          % in km
vJ   = sqrt(muS/aJ);            % Geschwindigkeit Jupiter heliozentr. km/s
RSOI = 48.2e06;                 % in km
RJ   = 69911;                   % Radius Jupiter in km
rPer = [10,15,20]*RJ;             % Perizentrumsabstand
gam1 = deg2rad(60);             % Anflugwinkel an SOI
v1   = 11.5;                    % Anfluggeschwindigkeit heliozentr. km/s

vinf = sqrt(v1^2+vJ^2-2*v1*vJ*cos(gam1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung ohne Boost 

exz1 = 1+rPer*vinf^2/muJ;
a1   = -rPer./sqrt((exz1-1).^2);
vPer = sqrt(muJ./rPer).*sqrt((exz1+1));
vPerk= sqrt(vinf^2+2*muJ./rPer);
thet1= 2*asin(1./(1+vinf^2*rPer/muJ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Berechnung mit Boost 

Delv = linspace(-1.5,4,100);       % Delta v in km/s

for k=1:length(Delv)
    for m=1:length(rPer)
        v2s(k,m)    = sqrt((vPer(m) + Delv(k))^2 -2*muJ/rPer(m));
        a2(k,m)     = -muJ/v2s(k,m)^2;
        exz2(k,m)   = (vPer(m)+Delv(k))^2*rPer(m)/muJ-1;
        thet2(k,m)  = pi-atan(sqrt(exz2(k,m)^2-1))-atan(sqrt(exz1(m)^2-1));
        Delth(k,m)  = rad2deg(thet2(k,m)-thet1(m));
        Delv2(k,m)  = v2s(k,m)-vinf;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphische Ausgabe I

titlestr = 'Gravity Assist mit Oberth Effekt am Jupiter ';
figure('name', titlestr)
grid on
subplot(1,2,1)
hold on
for m= 1:length(rPer)
    hp(m)= plot(Delv(:), Delv2(:,m),LC,Colors(m+1,:),LW,2);
    lgdstr(m,:) = sprintf('\\it r_{Per}/R_J\\rm = %4.1f',rPer(m)/RJ);
end
hp(m+1)=plot(Delv(:), Delv(:),LC,Colors(1,:),LS',':',LW,2);
lgdstr(m+1,:) = sprintf('\\Delta\\itv @ r = R_{SOI} ');
grid on
ttl = title('Änderung der Endgeschwindigkeit');
set(ttl,'FontSize',12, 'FontWeight','normal')
xlabel('\Delta \itv \rm in km/s','FontSize',12); 
ylabel('\Delta \itv_2'' \rm in km/s','FontSize',12);
hl1=legend(hp,lgdstr,...
                 'location','northwest');
legend box off
set(hl1,'FontSize',12)
set(gca,'FontSize',14)

subplot(1,2,2)
hold on
for m= 1:length(rPer)
    hp2(m)=plot(Delv(:), Delth(:,m),LC,Colors(m+1,:),LW,2);
end
grid on
ttl = title('Änderung der Ablenkwinkels');
set(ttl,'FontSize',12, 'FontWeight','normal')
xlabel('\Delta \itv \rm in km/s','FontSize',12); 
ylabel('Änderung der Ablenkwinkels in °','FontSize',12);
hl2=legend(hp2,lgdstr(1:3,:),...
                 'location','northeast');
legend box off
set(hl2,'FontSize',12)
set(gca,'FontSize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphische Ausgabe II


u=linspace(0,360,36000);

DelvSelect = [0,Delv(1),Delv(100)];
% Alles ink 
rR(1,:)  = a1(1)*(1-exz1(1)^2)./(1+exz1(1)*cosd(u(:)));
rR(2,:)  = a2(1,1)*(1-exz2(1,1)^2)./(1+exz2(1,1)*cosd(u));
rR(3,:)  = a2(100,1)*(1-exz2(100,1)^2)./(1+exz2(100,1)*cosd(u));
for m=1:3
 for k=1:length(u)
  xR(m,k) = rR(m,k)*cosd(u(k));
  yR(m,k) = rR(m,k)*sind(u(k));
  if xR(m,k)^2+yR(m,k)^2 > RSOI^2 || xR(m,k) > rPer(1) || (yR(m,k) < 0 && m>1)
      xR(m,k) = NaN;
      yR(m,k) = NaN;
  end
 end
end
xRSOI = RSOI*cosd(u);
yRSOI = RSOI*sind(u);
xJ    = RJ*cosd(u);
yJ    = RJ*sind(u);

% Bahnen

titlestr = 'Bahnen beim Gravity Assist mit Oberth Effekt';
figure('name', titlestr)
grid on
ttl = title(titlestr);
set(ttl,'FontSize',12, 'FontWeight','normal')
xlabel('\it x \rm in km','FontSize',12); 
ylabel('\it y \rm in km','FontSize',12);
hold on
axis equal

for m=1:3 
    hp3(m)=plot(xR(m,:),yR(m,:),'Color',Colors(m+1,:),'Linewidth',1,...
                'LineStyle',Style(1));
    lgdstr3(m,:) = sprintf('\\Delta\\it v\\rm = %4.1f km/s',DelvSelect(m));

end
hp3(4)=plot(xJ,yJ,'Color',Colors(5,:),LW,2,LS,Style(1));
hp3(5)=plot(xRSOI,yRSOI,'Color',Colors(5,:),LW,1,LS,Style(3));
grid on
lgdstr3(4,:) = ' Jupiter                  ';
lgdstr3(5,:) = ' SOI                      ';
hl2=legend(hp3,lgdstr3,'location','east');
legend box off
set(hl2,'FontSize',12)
set(gca,'FontSize',14)

% -------------------------------------------------------------------------
