% -------------------------------------------------------------------------
% KeplerVergleich.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programmerstellung  : Kenneth von Bünau
% -------------------------------------------------------------------------
% Lösung der Keplergleichung über verschiedene Iterationsverfahren 
% nach Kapitel 3.4.1 für verschiedene Exzentrizitäten 
% zwischen 0.1, 0.2,.....0.999999. 
% Vergleich mit Regula falsi und Newton-Raphson für verschiedene 
% Anfangswerte der mittleren Anomalie M.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
vals = 100; % Anzahl der betrachteten Punkte entlang der Bahn
M = linspace(0,2*pi,vals); % mittlere Anomalie

NrIter = 100; % Number of iterations  % Wähle bereich zwischen 5 und 1000
e = [0.2,0.5,0.8,0.9,0.95,0.999]; % Exzentrizitaeten

for n=1:length(e)
   C = [string(' \it{e}\rm=') num2str(e(n),'%5.3f')];
   lgdtitle(n,:) = strjoin(C);
end
Strich(1,:)= "-.";
Strich(2,:)= ":";
Strich(3,:)= "-";


% Iterationsverfahren------------------------------------------------------
E0 = ones(length(e),vals).*M;   % erster Wert
E0alt = sign(M-(M/2))*pi/2;     % alternative Anfangswerte

E = Iteration(M,e,E0,NrIter,vals);      % Standraditeration
Ealt = Iteration(M,e,E0alt,NrIter,vals);% Iterationen alternativer Startwert

% Regula falsi--------------------------------------------------------------
Ereg = RegulaFalsi(M,e,NrIter);

% Newton-Raphson------------------------------------------------------------
for m = 1:length(e)
    ENewton(m,:) = EAnom(M,e(m));
end

% Graphische Ausgabe--------------------------------------------------------

header=["Differenz E-M über Zeit für Iterationszahl = " , ...
        num2str(NrIter,'%3d')];
figure('Name',strjoin(header));
strNrIter = string(num2str(NrIter,'%3d'));
titles = ["Iterationsverfahren N="+strNrIter "Regula-falsi" ...
          "Newton-Raphson" "Iterationsverfahren N="+ ...
          strNrIter+newline+"Alternativer Startwert ±{\pi}/2"];
for k=1:4
    subplot(2,2,k);
    switch k
        case 1
            % plot Iterationsverfahren
            for m = 1:length(e)
                plot(M./(2*pi),E(m,:)-M,...
                'LineStyle', Strich(rem(m,3)+1), 'LineWidth',1);
                hold on;
            end
            % Beschriftung legende (Zuordnung der Exzentrititaeten)
            lgd=legend(lgdtitle,'location','southwest','numcolumns',2);
            legend box off;
        case 2
            % plot regula-falsi
            for m = 1:length(e)
                plot(M./(2*pi),Ereg(m,:)-M,...
                'LineStyle', Strich(rem(m,3)+1), 'LineWidth',1);
                hold on;
            end       
        case 3
           % plot Newton-Raphson
           for m = 1:length(e)
                plot(M./(2*pi),ENewton(m,:)-M,...
                'LineStyle', Strich(rem(m,3)+1), 'LineWidth',1);
                hold on;
            end
        case 4
            %plot Iteration alternativer Startwert
            for m = 1:length(e)
                plot(M./(2*pi),Ealt(m,:)-M,...
                'LineStyle', Strich(rem(m,3)+1), 'LineWidth',1);
                hold on;
            end
            lgd=legend(lgdtitle,'location','southwest','numcolumns',2);
    end
    ylim([-1 1]);
    ttl = title(titles(k));
    xlabel('Zeit/Umlaufperiode');
    ylabel('E-M');
    grid on;
    set(gca,'FontSize',14);
end
header='Lösung der Keplergleichung aus Iterationsverfahren';
text(10,11,header,'FontSize',18);

% Chart mit Abweichungen der Iteration zu Newton Raphson
header='Abweichung zu Newton-Raphson';
figure('Name',header);
ecc = [4 5 6];             % Ausgewählte Exzentrizitäten Laufindex m
iter = [3,5,25];           % Zahl von Iterationen Laufindex n
 
for m=1:3
    for n=1:3
        C = [lgdtitle(ecc(n)) " n="  num2str(iter(m),'%3d  ')]; 
        lgdtitle2(m,n,:) = strjoin(C);
    end
end

for k=1:2
subplot(1,2,k);
    switch k
        case 1
           %plot iterationsverfahren        
            for m=1:3
                for n = 1:3
                    E0 = ones(length(e),vals).*M;
                    E = Iteration(M,e(m),E0,iter(n),vals);
                    plot(M./(2*pi),E(ecc(m),:)-M-(ENewton(ecc(m),:)-M),...
                    'Color', Colors(m+1,:), 'LineStyle', Strich(n), ...
                    'LineWidth',1);
                    hold on;
                end
            end
            lgd=legend(lgdtitle2,'location','southeast','numcolumns',3);
         case 2
            %plot regula-falsi
            for m=1:3
                for n = 1:3
                    Ereg = RegulaFalsi(M,e,iter(n));
                    plot(M./(2*pi),Ereg(ecc(m),:)-M-(ENewton(ecc(m),:)-M),...
                    'Color', Colors(m+1,:), 'LineStyle', Strich(n),...
                    'LineWidth',1);
                    hold on;
                end
            end        
           lgd=legend(lgdtitle2,'location','southeast','numcolumns',3);
    end
    titles = ["Iterationsverfahren" "Regula-falsi" "Newton-Raphson" ...
    "Iterationsverfahren N="+strNrIter+newline+ ...
    "Alternativer Startwert ±pi/2"];
    ttl = title(titles(k));
    xlabel('Zeit/Umlaufperiode');
    ylabel('\Delta E');
    grid on;
    set(gca,'FontSize',16);
end

header='Lösung der Keplergleichung aus Iterationsverfahren';
text(10,11,header,'FontSize',18);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


% Funktionen --------------------------------------------------------------

function [E] = Iteration(M,e,E0,n,vals)
% berechnet nach Iterationsmethode mit n Iterationen und Startwert E0
% sowie vals M Werten
    E = E0;
    for m = 1:n % Iteration
        E = ones(length(e),vals).*M + (e.').*sin(E);
    end
end

function [Ereg] = RegulaFalsi(M,e,n)
% berechnet nach Regula falsi mit n Iterationen
Ereg = zeros(length(e),length(M));
for k = 1:length(e)
    for m = 1:length(M)
        % definiere anfängliches intervall [a,b] um M und M+e bzw M-e
        % sodass E auf jeden Fall enthalten ist
        if(M(m)>pi+0.5)
            a = M(m)-e(k);
            b = M(m)+2*pi/1000;
        elseif(M(m)<pi-0.5)
            a = M(m)-2*pi/1000;
            b = M(m)+e(k);
        else
            a = M(m)-e(k);
            b = M(m)+e(k);
        end
        
        % Iteration über umlauf mit n schritten
        for i = 1:n
            % interpoliere nach "regula falsi" einen wert linear im intervall
            c = (M(m)*(b-a)+(a*KeplerEq(b,e(k))...
                 -b*KeplerEq(a,e(k))))/(KeplerEq(b,e(k))-KeplerEq(a,e(k)));

            % probeweises errechnen des nächsten iterationsschritt mit
            % c als obere grenze
            d = (M(m)*(c-a)+(a*KeplerEq(c,e(k))...
                 -c*KeplerEq(a,e(k))))/(KeplerEq(c,e(k))-KeplerEq(a,e(k)));
            % setze neue intervallgrenzen
            % unterscheide dabei erste und zweite haelfte des umlaufs
            if((d<b && d>a && M(m)<pi))
                b = c;
            elseif(d<b && d>a)
                a = c;
            elseif(M(m)<pi)
                a = c;
            else
                b = c;
            end
        end
        
        Ereg(k,m) = (M(m)*(b-a)+(a*KeplerEq(b,e(k))-b...
            *KeplerEq(a,e(k))))/(KeplerEq(b,e(k))-KeplerEq(a,e(k)));
    end
end
end

function [M] =KeplerEq(E,e)
% Berechnet die Kepler-Gleichung
    M = E-e*sin(E);
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------