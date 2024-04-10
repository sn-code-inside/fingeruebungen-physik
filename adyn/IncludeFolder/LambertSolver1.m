% -------------------------------------------------------------------------
% LambertSolver1.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische FingerÃ¼bungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Parameter der Ellipsenbahn 
% für das Lambert Problem durch ein Bisektionsverfahren.
%
% -------------------------------------------------------------------------
% Inputs:         Description                     Range/Unit
%   p           - Vektor für p mit Startwerten         km
%   TOF         - TOF                                  s  
%   r1,r2       - Abstände P1, P2                      km 
%   gamma       - Winkel zwische Vektoren r1 und r2    rad
%   mu          - Gravitationsparameter                km^3/s^2 
%
% OutPuts:
%   p           - Vektor Iteration für p               km
%   ups1, ups2  - wahre Anomalien                      rad
%   pend        - Endwert p (semi-latus rectum)        km
%   eccend      - Elliptizitä e                        -
%   ups0        - Perihelwinkel                        rad
%   tFEnd       - Finaler Wert tF für Iteration        s
%
% -------------------------------------------------------------------------

% Lambert-Problem Iterationsverfahren Nr.1 
function [tFEnd,ups1,ups2,pend,eccend, ups0] = ...
          LambertSolver1(p,TOF,r1,r2,gamma,mu)
chi    = pi-gamma/2;  %Hilsfgröße in °
for k=1:2
    rho  = A(p(k),r1,r2,chi);
    sig  = B(p(k),r1,r2,chi);
    ecc(k) = sqrt(rho^2+sig^2);
    eta   = atan2(sig,rho);
    thet1 = chi-eta;
    thet2 = wrapTo2Pi(2*pi-(chi+eta));
    psi1  = 2*atan(sqrt((1-ecc(k))/(1+ecc(k)))*tan(thet1/2));
    psi2  = 2*atan(sqrt((1-ecc(k))/(1+ecc(k)))*tan(thet2/2));
    if psi2 <0 
        psi2 = psi2 +2*pi;
    end
    t1    = sqrt(p(k)^3/mu)*(psi1-ecc(k)*sin(psi1))/sqrt((1-ecc(k)^2)^3);
    t2    = sqrt(p(k)^3/mu)*(psi2-ecc(k)*sin(psi2))/sqrt((1-ecc(k)^2)^3);
    tF(k)   = t2-t1;
end
k=3;
while abs(tF(k-1)-TOF)/TOF > 0.001
    m(k) = (tF(k-2)-tF(k-1))/(p(k-2)-p(k-1)); 
    p(k) = (TOF - tF(k-2))./m(k) + p(k-2);
    rho  = A(p(k),r1,r2,chi);
    sig  = B(p(k),r1,r2,chi);
    ecc(k) = sqrt(rho^2+sig^2);
    eta    = atan(sig/rho);
    thet1 = chi-eta;
    thet2 = wrapTo2Pi(2*pi-(chi+eta));
    psi1  = 2*atan(sqrt((1-ecc(k))/(1+ecc(k)))*tan(thet1/2));
    psi2  = 2*atan(sqrt((1-ecc(k))/(1+ecc(k)))*tan(thet2/2));
    if psi2 <0 
        psi2 = psi2 +2*pi;
    end
    t1    = sqrt(p(k)^3/mu)*(psi1-ecc(k)*...
            sin(psi1))/sqrt((1-ecc(k)^2)^3);
    t2    = sqrt(p(k)^3/mu)*(psi2-ecc(k)*...
            sin(psi2))/sqrt((1-ecc(k)^2)^3);
    tF(k)   = t2-t1;
    tFEnd   = tF(k);
    ups1  = thet1;
    ups2  = thet2;
    pend   = p(k);
    eccend = ecc(k);
    ups0 = eta;
    k = k+1;
end 
end

% Bestimmung große Halbachse
function a = A(p,r1,r2,xi)
    a = (p/r1+p/r2-2)/2/cos(xi);
end

% Bestimmung kleine Halbachse
function b = B(p,r1,r2,xi)
    b = (p/r1-p/r2)/2/sin(xi);
end

