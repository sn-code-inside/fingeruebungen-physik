%%-------------------------------------------------------------------------
% AlternativeRaketenGl.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungend der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Berechnet die Raketengleichung nach einem diskreten Modell
%
% -------------------------------------------------------------------------

%% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter

mL   = 100;     % Masse Mensch (Nutzlast)
mB   = 100;     % Masse Boot (Strukturmasse)
mK   = 2.5;     % Masse Kokosnuss;
vK   = 10;      % Wurfgeschwindigkeit
Nmax = 100;     % Zahl der Kokosnüsse;
mT0  = Nmax*mK; % Gesamttreibstoffmasse


% Zeigt Annäherung an die Log-Funktionsberechnung
% mK fest, "Treibstoffmenge" mT = N*mK variiert
% mT0 2.5 ... 250 kg

for n=1:Nmax
    N(n) = n;
    V(n) = 0;
    mT(n) = n*mK;
    M0 = mL+mB+mT(n); 
    mu(n) = n*M0/mT(n);
    for k=1:n
        V(n) = V(n)+ vK*summand(k,mu(n));
    end
    Vlog(n) = vK*log(M0/(M0-mT(n)));
end

% Zeigt diskrete Annäherung an die Raketengleichung 
% mT0 fest, "Treibstoffmenge" pro "Zeitinterval" mT = MT0/N variiert
% mT 2.5 ... 250 kg

for n=1:Nmax
    N2(n) = n;
    mk(n)= mT0/n;
    V2(n) = 0;
    M0 = mL+mB+mT0; 
    mu2(n) = n*M0/mT0;
    for k=1:n
        V2(n) = V2(n)+ vK*summand(k,mu2(n));
        checksum(k) = 1/summand(k,mu2(n));
    end
    Vlog2(n) = vK*log(M0/(M0-mT0));
end



% Zeigt Annäherung an die Log-Funktionsberechnung
% mK fest, "Treibstoffmenge" mT = N*mK variiert
figure()
subplot(1,2,1)
plot(mT(:), V(:));
hold on
plot(mT(:), Vlog(:));
grid on
line([0 mT(length(N))],[vK vK])
ylabel('Geschwindigkeit Boot v_B','FontSize',14)
xlabel('"Treibstoffmenge" m_T(0) = N m_K','FontSize',14)
grid on;
h1 = title('Boot-Raketen-Analogie','FontSize',12);
set(h1,'FontSize',14,'FontWeight','normal'); 
h2=legend('Diskretes Model','Raketengleichung','v_K','location','southeast'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
axis([0, mT(length(N)), 0 1.2*vK]);
legend box off;
set(gca,'FontSize',16);

subplot(1,2,2)
plot(mk(:), V2(:));
hold on
plot(mk(:), Vlog2(:));
grid on
ylabel('Geschwindigkeit Boot v_B','FontSize',14)
xlabel('"Treibstoff pro Stufe" m_K = m_T(0) / N','FontSize',14)
grid on;
h1 = title('Boot-Raketen-Analogie','FontSize',12);
set(h1,'FontSize',14,'FontWeight','normal'); 
h2=legend('Diskretes Model','Raketengleichung','location','northwest'); 
set(h2,'FontSize',14,'FontWeight','normal'); 
legend box off;
set(gca,'FontSize',16);

%%-Ende Programm----------------------------------------------------------


%%-Funktionen-------------------------------------------------------------

function sumn = summand(n,mu)
  sumn = 1./(mu-n);
end

%%----------------------------------------------------------------------



