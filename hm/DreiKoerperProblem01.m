%%-------------------------------------------------------------------------
% DreiKoerproblem01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bewegung von drei Körpern, mit Anfansgbedingung auf 
% einem gleichseitigen Dreieck. 
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

RE  = 6371;         % in km

day  = 24 * 60 * 60; % Tag in s
year = 365.25 * day; % Jahr in s
AE   = 1.495978707e11;  % AE in km
G    = 6.6743e-11;      % G in m³ / (kg * s²)
m0   = 2e30;            % Massne in kg

% Simulationszeit und Zeitschrittweite [in s].
t_end = 50 * year;
dt    = 2 * day;


d = 1 * AE;     % Abstand des Körpers vom Schwerpunkt.
v0 = sqrt(G * m0 / d / sqrt(3)); %Geschwindigkeit von m0 (Betrag).

%  Anfangspositionen und -geschwindigkeitender Körper
Alpha = [0, 120, 240]+90;
vecr0 = d .* [cosd(Alpha); sind(Alpha)];
vecv0 = v0.* [-sind(Alpha); cosd(Alpha)];

figure()
hold on
for k=1:length(Alpha)
  plot(vecr0(1,k)/AE,vecr0(2,k)/AE,'o','markersize',8,'linewidth',8,...
       'color',Colors(k,:))
end
axis equal
axis([-2 2 -2 2]*3/4)
hold on
for k=1:length(Alpha)-1
    line([vecr0(1,1) vecr0(1,k+1)],[vecr0(2,1) vecr0(2,k+1)],'color','r')
end
line([vecr0(1,3) vecr0(1,2)],[vecr0(2,3) vecr0(2,2)],'color','r')
grid on

% %def dgl(t, u):
%     """Berechen Sie die rechte Seite der Differentialgleichung."""
%     r, v = np.split(u, 2)
%     r = r.reshape(n_körper, n_dim)
%     a = np.zeros((n_körper, n_dim))
%     für i im Bereich(n_körper):
%         für j im Bereich (i):
%             dr = r[j] - r[i]
%             gr = G / np.linalg.norm(dr) ** 3 * dr
%             a[i] += gr * m[j]
%             a[j] -= gr * m[i]
%     gibt np.concatenate([v, a.reshape(-1)]) zurück
% 
% 
% # Lege den Zustandsvektor zum Zeitpunkt t=0 fest.
% u0 = np.verketten((r0.reshape(-1), v0.reshape(-1)))
% 
% # Löse die Bewegungsgleichung numerisch.
% Ergebnis = scipy.integrate.solve_ivp(dgl, [0, t_max], u0, rtol=1e-9,
%                                    t_eval=np.arange(0, t_max, dt))
% t = Ergebnis.t
% r, v = np.split(Ergebnis.y, 2)
% 
% # Wandler und V in ein 3-dimensionales Array um:
% # 1. Index - Himmelskörper
% # 2. Index - Koordinatenrichtung
% # 3. Index - Zeitpunkt
% r = r.reshape(n_körper, n_dim, -1)
% v = v.reshape(n_körper, n_dim, -1)
% 
% # Berechne die verschiedenen Energiebeiträge.
% E_kin = 1/2 * m @ np.sum(v * v, Achse=1)
% E_pot = np.Null(t.Größe)
% für i im Bereich(n_körper):
%     für j im Bereich (i):
%         dr = np.linalg.norm(r[i] - r[j], Achse=0)
%         E_pot -= G * m[i] * m[j] / dr
% E = E_pot + E_kin
% dE_rel = (np.max(E) - np.min(E)) / E[0]
% print(f'Relative Energieänderung: {dE_rel:.2g}'
% 
% 
% # Erzeuge eine Figur und eine Achse für die Animation.
% fig = plt.figure()
% ax = fig.add_subplot(1, 1, 1)
% ax.set_xlabel('$x$ [AE]')
% ax.set_ylabel('$y$ [AE]')
% ax.set_xlim(-5, 5)
% ax.set_ylim(-5, 5)
% ax.set_aspect('gleich')
% ax.grid()
% 
% # Plotte für jeden Planeten die Bahnkurve.
% für ort, farbe in zip(r, farben):
%     ax.plot(ort[0] / AE, ort[1] / AE, '-',
%             Farbe=Farbe, Linienbreite=0,2)
% 
% # Erzeuge für jeden Himmelskörper einen Punktplot in der
% # entsprechenden Farbe und speichern Sie diese in der Liste.
% plots_himmelskoerper = []
% für Farbe in Farben:
%     plot, = ax.plot([], [], 'o', color=farbe)
%     plots_himmelskoerper.append(plot)
% 
% # Fügen Sie ein Textfeld für die Anzeige der verstrichenen Zeit hinzu.
% text_zeit = ax.text(-4.5, 4.5, '')
% 
% 
% def Aktualisierung(n):
%     "Aktualisiere die Grafik zum n-ten Zeitschritt."
%     für Plot, Ort in zip(plots_himmelskoerper, r):
%         plot.set_data(ort[:, n].reshape(-1, 1) / AE)
%     text_zeit.set_text(f'{t[n] / jahr:.2f} Jahre')
%     return plots_himmelskoerper + [text_zeit]
% 
% 
% # Erzeuge das Animationsobjekt und starte die Animation.
% ani = mpl.animation.FuncAnimation(Abb., aktualisieren, Frames=t.size,
%                                   Intervall=30, blit=True)
% 
% %% VorwÃ¤rtsrechnung Orbit aus Bahnparametern
% 
% % Berechnung mit SSO Parametern
% eP     = 0.00;
% OmegaP = deg2rad(100.1213);
% omegaP = deg2rad(0);
% MP     = deg2rad(wrapTo360(0));  
% K      =  3;  %Zyklus in Tage
% N      = 43;  %Zyklus in Umrundungen der Erde
% D_La1  = -360*K/N;  %geforderte Verschiebung aufsteigender Knoten
% %Berechnet aus Bedingung N*(OmegaDot - ThetaDot)+TP = -K*360°
% aP     = 7153.123;   %in km
% J2     = 1.083e-3;  
% iP     = deg2rad(98.498);
% NPoints= 201;
% 
% %rate change of Omega per second 
% OmegaDot  = -1.5*sqrt(GME/aP^3)*J2*cos(iP)*RE^2/aP^2;  % per second
% OmegaDotd = OmegaDot*24*60*60;   %rate change of Omega per day 
% 
% %rate change of omega per day 
% omegaDot  = -0.75*sqrt(GME/aP^3)*(1-5*(cos(iP))^2)*J2*RE^2/aP^2; % per second 
% omegaDotd = omegaDot*24*60*60;   %rate change of Omega per day 
% 
% %rate change of mean anomaly
% nDot  = -0.75*sqrt(GME/aP^3)*(1-3*(cos(iP))^2)*J2*RE^2/aP^2; % per second
% nDotd = nDot*24*60*60;    %rate change of Omega per day
% n0d   = sqrt(GME/aP^3)*86400;
% 
% %rate change of Theta
% ThetaDotd = 360.9856473*2*pi/360; %rad per day
% ThetaDot  = ThetaDotd/24/60/60;   %rad per second
% MP        = wrapTo360(0);  % mean anomaly 
% 
% %Umlaufzeiten
% TPn    = 2*pi*sqrt(aP^3/GME)/86400;  %Umlaufzeit Kepler-Bahn in d
% TPkor  = 2*pi/(n0d+nDotd+omegaDotd); %geforderte Umlaufzeit für SSO
% 
% %Check for closure after N=3 and K=43
% % x1 = N*(OmegaDotd-ThetaDotd)*2*pi/(n0d+nDotd+omegaDotd);
% % x2 = -K*2*pi;
% % Parameterübergabe
% 
% TSpan = 6/24; %Zeitspanne in h
% 
% % Berechnung Bahndaten zeitpunkt T1
% dt1 = datetime('2020-01-01 00:00:00');
% T1  = juliandate(dt1) - TSpan/24; % Julianisches Datum  ET
% T1E = T1 + TSpan/24 + TPkor;
% T_vector1 = linspace(T1,T1E,NPoints);  %in julianischen Tagen udn Bruchteilen
% SatSSO.BaPa =[aP , eP, iP , MP, OmegaP, omegaP];
% SatSSO.Name = 'SSO';
% SatData1  = SatPQR_perturbed_orbit(T_vector1, SatSSO, GME, OmegaDotd, omegaDotd, nDotd);
% 
% % Berechnung Bahndaten zeitpunkt T2
% T2  = T1 + 42*TPkor  - TSpan/24;  % Julianisches Datum  ET
% T2E = T1 + 43*TPkor  + TSpan/24;
% MP = MP + (n0d + nDotd + omegaDotd)*(T2-T1) ;
% OmegaP = OmegaP+OmegaDotd*(T2-T1);
% T_vector2 = linspace(T2,T2E,NPoints);  %in julianischen Tagen udn Bruchteilen
% SatSSO.BaPa =[aP , eP, iP , MP, OmegaP, omegaP];
% SatData2  = SatPQR_perturbed_orbit(T_vector2, SatSSO, GME, OmegaDotd, omegaDotd, nDotd);
% 
% % Umrechnung in geographische Breite und Länge
% lat1(:) = rad2deg(SatData1.el);
% theta1  = GMSTsat(T_vector1);
% lon1(:) = wrapTo180(rad2deg(SatData1.az)-theta1);
% lat2(:) = rad2deg(SatData2.el);
% theta2  = GMSTsat(T_vector2);
% lon2(:) = wrapTo180(rad2deg(SatData2.az)-theta2);
% 
% 
% %% Graphik Ground Track
% 
% titlestr = strcat(SatSSO.Name,' Ground Track');
% figure('name',titlestr);
% gx = geoaxes;
% hp(1)=geoplot(gx, lat1(:),lon1(:),'d',...
%         'MarkerSize',2,'Color', Colors(2,:),'LineWidth', 3);
% hold on
% hp(2)=geoplot(gx, lat2(:),lon2(:),'+',...
%         'MarkerSize',2,'Color', Colors(4,:),'LineWidth', 2);
% geobasemap(gx,'bluegreen')
% geolimits('manual') 
% geolimits([-90 90],[-180 +180])
% hp1 = title(titlestr,'FontSize',12);
% legend(hp,'1. Orbit', '42. und 43. Orbit', 'location', 'northeast')
% % legend box off
% set(hp1,'FontSize',14,'FontWeight','normal'); 
% set(gca,'FontSize',14);
% 
% %% PrintOut
% 
% 
% fprintf('\n');
% fprintf('Bahnparameter fÃ¼r Bahnberechnung nach Kepler\n');
% fprintf('\n Zeit  = %s', dt1);
% fprintf('\n');
% fprintf('\n a     = %8.2f km  (GroÃe Halbachse)', aP);
% fprintf('\n h     = %8.2f km  (mittlere Höhe)', aP-RE);
% fprintf('\n e     = %8.5f     (ExzentrizitÃ¤t)', eP);
% fprintf('\n i     = %8.2f     (Inklination)', iP);
% fprintf('\n Omega = %8.2f     (Rektasz. aufst. Knoten)', rad2deg(OmegaP));
% fprintf('\n omega = %8.2f     (Argument PerigÃ¤um)', rad2deg(omegaP));
% fprintf('\n M     = %8.2f     (Mittlere Anomalie)', rad2deg(MP));
% fprintf('\n TPn   = %8.2f min (Umlaufzeit Kepler-Bahn)', TPn*24*60);
% fprintf('\n TPkor = %8.2f min (Umlaufzeit Sonnensynchrone Bahn)', TPkor*24*60);
% fprintf('\n');
% fprintf('\n D_La1 =%+8.2f min (Verschiebung aufsteigender Knoten gefordert)', D_La1);
% fprintf('\n');
% fprintf('\n');

%% Ende Programm


