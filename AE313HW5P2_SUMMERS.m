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



