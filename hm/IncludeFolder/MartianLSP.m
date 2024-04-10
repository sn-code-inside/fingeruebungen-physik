% -------------------------------------------------------------------------
% MartianLSP.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Derivation of LSP from ecliptic position of Pole Vector
% -------------------------------------------------------------------------

clc
close all 
Colors=GetColorLines;

% Initialisierung 
dt1 = datetime('2000-01-02 12:00:00');
% Julianische Startdatum in UT
T1      = juliandate(dt1);
epsE    = EpsErde(T1);  
epsErad = deg2rad(epsE);  

Aequi = 'Datum';
Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T1,Aequi);

% position of planet polar axis in geocentric equatorial coords [3]
alphaA(4)  = 317.681;
deltaA(4)  =  52.886;
alphaA(3)  =   0.0;
deltaA(3)  =  90.0;
fprintf('\n');

for PlNr=3:4  
    
    % orbit parameters
    iP     = BaPa.iP(PlNr);      % orbital inclination mars
    OmegaP = BaPa.OmegaP(PlNr);  % longitude of node
    PiqP   = BaPa.OmegaP(PlNr);  % longitude of perihelion
    omega  = PiqP-OmegaP;        % argument of perihelion
    M0     = BaPa.M0P(PlNr);     % mean anomalie at epoch
    alpha0 = alphaA(PlNr);
    delta0 = deltaA(PlNr);
 
    % convert to geocentric ecliptical coords
    axis_b = asind(-cosd(delta0)*sind(alpha0)*sind(epsE)+ ...
         sind(delta0)*cosd(epsE));
    axis_l = wrapTo360(atan2d(cosd(delta0)*sind(alpha0)*cosd(epsE)+ ... 
                sind(delta0)*sind(epsE),cosd(delta0)*cosd(alpha0)));    
    % rotation matrices
    Rz_O  = [cosd(OmegaP) -sind(OmegaP) 0; sind(OmegaP) cosd(OmegaP) 0;...
            0 0 1];
    Rz_mO = [cosd(-OmegaP) -sind(-OmegaP) 0; ...
             sind(-OmegaP) cosd(-OmegaP) 0;0 0 1];
    Rz_o  = [cosd(omega) -sind(omega) 0; sind(omega) cosd(omega) 0;0 0 1];
    Rx_i  = [1 0 0;0 cosd(iP) -sind(iP); 0 sind(iP) cosd(iP)];
    Rx_mi = [1 0 0;0 cosd(-iP) -sind(-iP); 0 sind(-iP) cosd(-iP)];
     %convert polar axis to cartesian coords
    p                = [0,0,0];
    [p(1),p(2),p(3)] = sph2cart(axis_l,axis_b,1);
    % calculate angular momentum vector of orbit in ecliptical coordinates
    L_bahn = [0,0,1];        % in cartesian orbit coordinates
    L = Rz_mO*Rx_i*L_bahn';  % in ecliptical cartesian coordinates
    % calculate the vernal direction in the ecliptic
    vernal_ecliptic = cross(p',L); 
    % calculate the vernal direction with respect to perihelion of mars
    vernal = Rz_o*Rx_mi*Rz_O*vernal_ecliptic; % in cartesian ecliptic 
    % convert to spherical coordinates and display longitude LSP
    [LSP,b,c]=cart2sph(vernal(1),vernal(2),vernal(3));
    LSP = wrapTo360(rad2deg(LSP));
    fprintf('%s \t alpha0: % 8.3f°\t delta0  : % 8.3f° \t axis_l: % 8.3f°\t axis_b: % 8.3f° \n', ...
            string(BaPa.Name(PlNr)), alpha0 , delta0,...
            axis_l, axis_b);
    fprintf('%s \t LSP   : % 8.3f°\t ALphaFMS: % 8.3f° \n', ...
            string(BaPa.Name(PlNr)), LSP , wrapTo360(LSP+M0));
end
fprintf('\n');

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------