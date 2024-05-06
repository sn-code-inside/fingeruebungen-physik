% -------------------------------------------------------------------------
% RutherfordAg.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Rutherford-Streung Ag auf Au
%
% Programm berechnet Wirkungsquerschnitte für
% Rutherford-Streuung mit und ohne Rückstoß
%
% Streuteilchen: Silber
% Streuzentrum:  Gold
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];
Marker = ['o','d','o','s','+'];

% Variablen, Konstanten 
u     = 1.660539e-27;       % Atomare Masseneinheit kg
eq    = 1.602e-19;          % Elementarladung in As
eps0  = 8.854e-12;          % Dielektrizitaetskonstante As/Vm
aB    = 5.29177210e-11;     % Bohrradius in m
RZ = 0.0;                   % Zenmtrumsposition
ZT = 47; ZZ = 79;           % zT=Teilchen, zZ=Streuzentrum Ladungen
kappa  = eq^2*ZZ*ZT/4/pi/eps0;
mT  = 109*u;                % Teilchenmasse in atomic units
mZ  = 197*u;                % Streuzentrumsmasse in atomic units
En_eV   = [2e4,1e4,5e4];    % Energie in eV 
En  = En_eV/6.242e+18;      % Umrechnung in Ws
vinf   = sqrt(2*En/mT);     % Initialgeschwindigkeit im Unendlichen in m/s
L=25;                       % Betrachtungsfenster in units of b
b   = 3*aB;                 % Stoßparameter
x10 = 50*aB;                % Initiale x Position --> unendlich
tmax = 2*x10/vinf(2);
y10 = b;                    % Initiale y Position = b und Geschw.
d1x0  =-sign(x10)*sqrt(2/mT*(En-kappa/sqrt(x10^2+y10^2))); 
d1y0=0.0;                   % Null y-speed 
for k = 1:3                 % AB
    AB(k,:) =[x10,d1x0(k),y10,d1y0,0,0,0,0];
end
tv = linspace(0,tmax,100);

% Berechnung 
for k=1:3
    [tout,Yout]   = asteroid_path(L, kappa, mT, mZ, tv, AB(k,:));
    Data(k).xT = Yout(:,1)/aB;
    Data(k).yT = Yout(:,3)/aB;
    Data(k).xZ = Yout(:,5)/aB;
    Data(k).yZ = Yout(:,7)/aB;
    tend(k) = tout(end)/1e-15;
end
lx1 = length(Data(1).xT);


%% 
% Graphische Ausgabe

figure('name','Rutherford-Streuung Ag an Au' )
axis ([-L L -L L])
grid on
ttl = title('Rutherford-Streuung Ag an Au' );
set(ttl,'FontSize',14, 'FontWeight','normal')
xlabel('\it x \rm in \it a_B','FontSize',13); 
ylabel('\it y \rm in \it a_B','FontSize',14);
hold on
for ki=1:1:lx1 % all points animated
    p(1)=plot(0,0,'s','MarkerFaceColor', 'w','LineWidth',2,...
             'MarkerEdgeColor',Colors(4,:),'MarkerSize',10);%Zentrum
    plot(x10/aB ,b/aB ,'o','MarkerFaceColor',Colors(4,:),...
        'MarkerEdgeColor',Colors(4,:),'MarkerSize',6); %Teilchen
    p(2)=plot(L ,b/aB ,'o','MarkerFaceColor','w','LineWidth',2,...
        'MarkerEdgeColor',Colors(4,:),'MarkerSize',8); %Teilchen
    p(3) = plot(Data(1).xT ,Data(1).yT ,'Color', Colors(4,:),...
         'LineWidth',2,'LineStyle',Style(1)); %plot yT versus xT
    plot(Data(1).xZ ,Data(1).yZ ,'Color', Colors(4,:),...
          'LineWidth',2,'LineStyle',Style(2)); %plot yZ versus xZ
    plot(Data(1).xT(ki) ,Data(1).yT(ki) ,'o', 'MarkerFaceColor', Colors(4,:),...
         'MarkerEdgeColor',Colors(4,:),'MarkerSize',6); %Teilchen
    plot(Data(1).xZ(ki) ,Data(1).yZ(ki) ,'s','MarkerFaceColor', Colors(4,:),...
         'MarkerEdgeColor',Colors(4,:),'MarkerSize',7); %Zentrum
    pause(0.2)
end
for k =2:3
    lx = length(Data(k).xT);
    ptemp = plot(Data(k).xT ,Data(k).yT ,'Color', Colors(k,:),...
     'LineWidth',1,'LineStyle',Style(1)); % plot yT versus xT
    if k == 2 
    p(4) = ptemp; end
    ptemp = plot(Data(k).xZ ,Data(k).yZ ,'Color', Colors(k,:),...
      'LineWidth',1,'LineStyle',Style(2)); % plot yZ versus xZ
    if k == 3 
    p(5) = ptemp; end
    plot(Data(k).xT(lx) ,Data(k).yT(lx) ,'o','MarkerFaceColor',Colors(k,:),...
         'MarkerEdgeColor',Colors(k,:),'MarkerSize',6); % Teilchen
    plot(Data(k).xZ(lx) ,Data(k).yZ(lx) ,'s','MarkerFaceColor',Colors(k,:),...
         'MarkerEdgeColor',Colors(k,:),'MarkerSize',7); % Zentrum
end

% Eingangsasymptote
line([-L,L],[b/aB,b/aB],'Color',Colors(9,:),'LineStyle',Style(3),...
     'LineWidth',1); %Eingangsasymptote
sa=(Data(1).yT(lx1)-Data(1).yT(lx1-1))/(Data(1).xT(lx1)-Data(1).xT(lx1-1)); 
% Ausgangsasymptote  
xa=-L/3:L;
ya=Data(1).yT(lx1) +sa*(xa-Data(1).xT(lx1) );
plot(xa,ya,'Color',Colors(9,:),'LineStyle',Style(3),...
     'LineWidth',1);                
% Berechne r(t) and rmin
r=sqrt((Data(1).xT-Data(1).xZ).^2+(Data(1).yT-Data(1).yZ).^2);
[rmin, imin] = min(r);
line([Data(1).xT(imin) ,Data(1).xZ(imin) ],...
     [Data(1).yT(imin) ,Data(1).yZ(imin) ],...
     'Color',Colors(9,:),'LineStyle',Style(3),'LineWidth',1);
line([-L/3,L],[0,0],'Color',Colors(15,:),'LineStyle',Style(3),...
     'LineWidth',1);
line([0,0],[-L,L/3],'Color',Colors(15,:),'LineStyle',Style(3),...
     'LineWidth',1);


 
% Parameterausgabe
str1(1,:) = string(cat(2,' \it b   \rm  = ',...
            num2str(b/aB ,'%4.2f'),' a_B  '));
str1(2,:) = string(cat(2,' \it r_{min}\rm = ',...
            num2str(rmin ,'%4.2f'), ' a_B  '));
for k = 1:2 
text(-L*(1-0.2),-L+k*3,str1(k,:),'FontSize',14,'FontWeight','normal');
end

% Legende
for k=1:3
    strEn(k,:) = string(cat(2,...
         ' \it E_{inf} \rm = ',num2str(En_eV(k)/1000,'%3.1f'),' keV',...
         '  t_{end} = ',num2str(tend(k) ,'%4.2f'),' fs  '));
end
h1=legend(p(1:5),'Au - Zentrum','Ag - Teilchen',strEn(1),strEn(2),...
                 strEn(3),'location','bestoutside');
legend box off
axis square

ttl = title('Rutherford-Streuung Ag an Au' );
set(ttl,'FontSize',14, 'FontWeight','normal')
set(h1,'FontSize',14)
set(gca,'FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
% Funktionen

function [tout,Yout] = asteroid_path(L, kappa, mT, mZ, tv, AB)
    opt = odeset('AbsTol',1.e-9,'RelTol',1.e-7,'Events',@events);
    [tout,Yout,te,ye,ie]=ode45(@(t,Y)dgl_ruther(t,Y,kappa,mT,mZ),tv,AB,opt);
    function [value,isterminal,direction] = events(t,Y)
        aB    = 5.29177210e-11;     % Bohrradius in m
        r1=sqrt(Y(1).^2 + Y(3).^2);
        value = (r1/aB - L);        % detect distance 
        isterminal = 1;             % stop the integration
        direction  = 1;             % negative direction
    end
end


% DGL-System
function dY = dgl_ruther(t, Y, kappa, mT, mZ)
    % Y(1,2,...8)=x1,v1x,y1,v1y,x2,v2x,y2,v2y
    % es muss ein Spaltenvektor zurückgegeben werden 
    dY     =  zeros(8,1); 
    r      = sqrt((Y(1)-Y(5)).^2 + (Y(3)-Y(7)).^2);
    dY(1)  =  Y(2);
    dY(2)  = kappa*(Y(1)-Y(5))./r.^3/mT;
    dY(3)  =  Y(4);
    dY(4)  = kappa*(Y(3)-Y(7))./r.^3/mT;
    dY(5)  =  Y(6);
    dY(6)  = -kappa*(Y(1)-Y(5))./r.^3/mZ;
    dY(7)  =  Y(8);
    dY(8)  = -kappa*(Y(3)-Y(7))./r.^3/mZ;
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
