% Name: Shepard Summers
% Date: 10/30/2024
% Professor: Hao Peng
%
% Program purpose: To use Radar method given angles
%
% Assumptions: Coplanar 

clear; clc;

m = 398600; % km^3/s^
rE = 6378; % km
omE = [0, 0, 7.2921 * 10^-5].'; % rad/s

% Input variables

% Azimuth and Eleveation deg --> rad
A = 36 * (pi / 180); 
ADot = 0.590 * (pi / 180);
a = 36.6 * (pi / 180);
aDot = -0.263 * (pi / 180);

% Range 
p = 988; 
pDot = 4.86;

% Radar station
theta = 40 * (pi/180);
lat = 35 * (pi/180);
alt = 1;
rR = rE + alt;

% p hat in ijk frame
pHat_ijk = [cos(a)*sin(A), cos(a)*cos(A), sin(a)].';

% convert p hat to IJK
psi = [-sin(theta), -sin(lat)*cos(theta), cos(theta)*cos(lat);
    cos(theta), -sin(lat)*sin(theta), sin(theta)*cos(lat);
    0, cos(lat), sin(lat)];

pHat_IJK = psi * pHat_ijk;

% calculate r
rEVect_IJK = [rE*cos(lat)*cos(theta), rE*cos(lat)*sin(theta), rE*sin(lat)].';

pVec = p * pHat_IJK;

r = rEVect_IJK + pVec;

% calculate v
pDotHatRel_ijk = [-sin(a)*aDot*sin(A) + cos(a)*cos(A)*ADot, -sin(a)*aDot*cos(A) - cos(a)*sin(A)*ADot, cos(a)*aDot].';

pDotHatRel_IJK = psi * pDotHatRel_ijk; 

pDotHat_IJK = pDotHatRel_IJK + cross(omE, pHat_IJK);

rDot = cross(omE, r);

v = rDot + pDot*pHat_IJK + p*pDotHat_IJK;

fprintf("r");
disp(r);
fprintf("v");
disp(v);

% Problem 3 =================================

rMag = sqrt(sum(r.^2));

% Calculate h
hVec = cross(r, v);
h = sqrt(sum(hVec.^2));

% Caclulate e
eVec = cross(v, hVec)/m - r/rMag;
e = sqrt(sum(eVec.^2));

% Calculate mean anomoly
maRad = acos((h^2)/(m*rMag*e) - 1/e);
if dot(r, v) < 0  
    % angle correction due to unknown quadrent
    maRad = 2 * pi - maRad;
end
ma = maRad * 180/pi;

% Node line needed for angles
nodeLine = cross([0,0,1], hVec);
nodeMag = sqrt(sum(nodeLine.^2));
    
% Right Ascension of the Ascending Node
ra = acos(nodeLine(1) / nodeMag); 
if nodeLine(2) < 0
    % angle correction due to unknown quadrent
    ra = 2 * pi - ra;
end
ra = (180 / pi) * ra;

% Inclination
inc = acos(dot(hVec, [0, 0, 1]) / h);
inc = (180 / pi) * inc;
    
% Argument of Perigee
ap = acos(dot(nodeLine, eVec) / (nodeMag * e));
if eVec(3) < 0
    % angle correction due to unknown quadrent
    ap = 2 * pi - ap;
end
ap = (180 / pi) * ap;

a = ((h^2) / m) / (1-e^2);

% Output data
fprintf("h");
disp(h);
fprintf("e");
disp(e);
fprintf("theta");
disp(ma);
fprintf("RAAN");
disp(ra);
fprintf("INC");
disp(inc);
fprintf("AP");
disp(ap);
fprintf("a");
disp(a);



