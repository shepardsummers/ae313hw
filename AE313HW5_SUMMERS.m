% Name: Shepard Summers
% Date: 10/30/2024
% Professor: Hao Peng
%
% Program purpose: To use Gibbs method for 3 given position vectors
%
% Assumptions: Coplanar 

clear; clc;

m = 398600; % Mew for Earth in km^3/s^2

% Position vectors
r1 = [5887, -3520, -1204];
r2 = [5572, -3457, -2376];
r3 = [5088, (-3289+10), -3480];
r = {r1, r2, r3};

% Position magnitudes 
r1Mag = sqrt(sum(r1.^2));
r2Mag = sqrt(sum(r2.^2));
r3Mag = sqrt(sum(r3.^2));
rMag = [r1Mag, r2Mag, r3Mag];

% N, D, S vectors
N = r1Mag * cross(r2, r3) + r2Mag * cross(r3, r1) + r3Mag * cross(r1, r2);
D = cross(r1, r2) + cross(r2, r3) + cross(r3, r1);
S = (r2Mag - r3Mag) * r1 + (r3Mag - r1Mag) * r2 + (r1Mag - r2Mag) * r3;

% N, D magnitudes
nMag = sqrt(sum(N.^2));
dMag = sqrt(sum(D.^2));

v = cell(1);
hMag = zeros(1);
eMag = zeros(1);

for i = 1:3
    vI = sqrt(m / (nMag * dMag)) * (cross(D, r{i})./rMag(i) + S);
    v{i} = vI;

    hI = cross(r{i}, v{i});
    hMag(i) = sqrt(sum(hI.^2));

    eI = (cross(r{i}, v{i}) / m) - (r{i} / rMag(i));
    eMag(i) = sqrt(sum(eI.^2));
    
    
end









fprintf("+-+-+ OUTPUTS +-+-+")
fprintf("\n Velocity vectors: ")
disp(v);
fprintf("\n h magnitude: ")
disp(hMag);
fprintf("\n e magnitude: ")
disp(eMag);

