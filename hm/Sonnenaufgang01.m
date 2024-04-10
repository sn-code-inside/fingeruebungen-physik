% -------------------------------------------------------------------------
% Sonnenaufgang01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die Aufgangs- und Untergangszeiten der Sonne
% für verschiedene geographische Orte (Breite < 60°) und die Taglängen
% iterativ aus der Keplergleichung für das ganze Jahr und in der Nähe des 
% Jahreswechsels bzw. Mitsommers. Es wird sichtbar, dass der kürzeste Tag  
% des Jahres nicht mit dem  spätesten Sonnenaufgang 
% und frühesten Sonnenuntergang zusammenfällt.
%
% Achten Sie bei der Ausführung des Programms auf Ihre Kommandozeile
% (Command Window) und bestätigen Sie mit "J", falls die Daten für
% Oberkochen verwendet werden sollen. Möchten Sie dieses Programm für einen
% anderen Ort verwenden, wählen Sie "N".
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

pi2 = 2*pi;                  
lambdaOko=10.0990193; %Geographische Koordinaten für Oberkochen
phiOko=48.7854368;
hr=deg2rad(-50/60); %Höhenkorrektur(Parallaxe und scheinbarem Durchmesser)
pi2 = 2*pi;

%_________________________________________________________________________
% Initialisierung

% Abfrage der Eingabedaten

Oko     = input('\Oberkochen J/N  :                 ... ','s');
if Oko == 'J' || Oko =='j' 
    Oko='j';
    lambda=lambdaOko;
    phi=deg2rad(phiOko);
    Zone=1; 
else
    Oko='n';
    phi     = input('\Geographische Breite +-BB.bb      ... ');
    lambda  = input('\Geographische Laenge +-LL.ll      ... ');
    Zone    = input('\Zeitzone inkl. DST +-UT           ... ');
    phi=deg2rad(phi);
    sinphi=sin(phi);
    cosphi=cos(phi);
end
MOZ= lambda/15;  % mittlere Ortszeit (s. Kap. Raum und Zeit)
eps0=deg2rad(23.43929111);      % Schiefe der Ekliptik

dt1 = datetime('2000-01-01 12:00:00');
T = juliandate(dt1); % Julianisches Datum
MJuDa = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
 
GMThour = 12 - Zone;

zeitdiff=zeros(1,400);
Deklination=zeros(1,400);
Alt=zeros(1,400);
tauw=zeros(1,400);  
Timey=zeros(1,400);
kulmi=zeros(1,400);
taglaenge=zeros(1,400);
auf=zeros(1,400);
unter=zeros(1,400);

%------------------------------------------------------------------------------

%Beginn Rechnung

T_v=T:(T+399);  %Bewusst länger als ein Jahr, da Berechnung über 
                     %Jahreswechsel
MJuDa_v = MJuDa:(MJuDa+399);
Timey = 1:400;

% Berechnung Rektaszension und Deklination über Mittelpunktsgleichung;
[RA_v,Dec_v]=KeplerSonne(T_v,eps0); 
% Stundenwinkel in Bogenmaß
tau1=(GMST(MJuDa_v)+deg2rad(lambda)-RA_v);  
% Stundenwinkel im Zeitmaß
tauw = rad2deg(tau1)/15 +(12-GMThour)-MOZ; % in Stunden
tauw = 24*wrapToPi(pi2*tauw/24)/pi2;      % in Stundem [0..24]
q(:)=(sin(hr)-sin(phi)*sin(Dec_v))./(cos(phi)*cos(Dec_v));
kulmi=12-tauw-MOZ+Zone;  %Kulmination

for k=1:400  
    if abs(q(k)) < 1
        zeitdiff(k)    = acosd(q(k))/15; 
        auf(k)    = kulmi(k) - zeitdiff(k);
        unter(k)  = kulmi(k) + zeitdiff(k);
        taglaenge(k)   = 2*zeitdiff(k);
    else
        if q(k)<0 
            taglaenge(k)= 23.95;
         else
            taglaenge(k)= 0.05; 
            kulmi(k)=NaN;  %Kulmination
        end
        zeitdiff(k) = acosd(q(k))/15; 
        auf(k)    = NaN; 
        unter(k)  = NaN; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphische Ausgabe
phistr=sprintf(' %+4.1f°',rad2deg(phi));
lstr=sprintf(' %+4.1f°',lambda);    

% Bild 1
% SA, SU und Taglänge im Verlauf eines Jahres
% Polarsommer und Polarwinter werden ausgeblockt
figure();
% subplot(1,2,1) 
for k=1:12
    XDataTick(k) = datetime(2020,k,1);
end
XDataTick(13) = datetime(2021,1,1);
xData=datetime('2020-01-01') + caldays(1:400);
plot(xData, auf,'Color',Colors(2,:),'LineWidth',2); 
hold on;
plot(xData, kulmi,'Color',Colors(5,:),'LineWidth',2);
plot(xData, unter,'Color',Colors(3,:),'LineWidth',2);
plot(xData, taglaenge,'Color',Colors(2,:),'LineWidth',2,'LineStyle',':');
if Oko == 'J' || Oko =='j'  
    txt = ['SA,SU,Taglänge in Oberkochen'];
else
    txt = ['SA,SU,Taglänge'];
end
text(5,22.66, txt,'FontSize',18);
txt= ['Breite ',phistr,' und Länge  ',lstr];
text(5,21.5,txt,'FontSize',18);
grid on;
ax = gca;
ax.XTick = XDataTick;
datetick('x','mmm','keepticks')
lgd=xlabel('Monat');
lgd.FontSize=18;
ylim([0,24]);
lgd=ylabel('Zeit in MOZ bzw. Länge in h');
lgd.FontSize=18;
yticks manual;
yticks([0 3 6 9 12 15 18 21 24]);
yticklabels({'00:00' '03:00' '06:00' '09:00' '12:00' '15:00' '18:00' '21:00' '24:00'});
lgd=legend('SA','Mittag','SU','Taglänge','location','southeast');
legend boxoff;
lgd.FontSize=18;
hold on;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

% Bild 2
% Darstellung der Situation um die Wintersonnenwende und Periheldurchgang
if (rad2deg(phi) < (90-rad2deg(eps0))) && (phi > 0) 
    figure()
    % subplot(1,2,2);
    for k=1:10
        XDataTick2(k) = datetime(2020,11,1+10*(k-1));
    end
    xData=datetime('2020-01-01') + caldays(1:400);
    yyaxis left
    plot(xData, taglaenge,'Color',Colors(2,:),'LineStyle',':','LineWidth',2);
    hold on;
    plot(xData, auf,'Color',Colors(2,:),'LineStyle','-','LineWidth',2);
    hold on;
    ax=gca;
    ylim([auf(365)-2,auf(365)+2]);
    xlim([xData(332) xData(400)]);
    ax.YColor = Colors(2,:);
    Colors(1,:);
    lgd=ylabel('SA-Zeit in MOZ bzw. Taglänge in h');
    lgd.FontSize=18;
    lgd.Color=Colors(2,:);
    ylim([7 10]);
    yticks manual;
    yticks([7 7.5 8 8.5 9 9.5 10]);
    yticklabels({'7:00','7:30','8:00','8:30','9:00','9:30','10:00'})
    yyaxis right
    plot(xData, unter,'Color',Colors(3,:),'LineWidth',2);
    mid=round(unter(365));
    ylim([mid-2,mid+2]);
    lgd=ylabel('SU-Zeit in MOZ');
    lgd.FontSize=18;
    lgd.Color=Colors(3,:);
    ax.YColor = Colors(3,:);
    lgd=legend('Taglänge','SA','SU','location','southeast');
    legend boxoff;
    lgd.FontSize=18;
    txt= 'SA/SU und Tageslänge am Jahreswechsel';
    text(325,18.75,txt,'FontSize',18);
    grid on;
    ax = gca;
    for k=1:6
        XDataTick3(k) = datetime(2020,11,1+20*(k-1));
    end
    ax.XTick = XDataTick3;
    datetick('x','dd.mm.','keepticks');
    lgd=xlabel('Tag');
    lgd.FontSize=18;
    set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------