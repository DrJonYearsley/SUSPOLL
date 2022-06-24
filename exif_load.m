% Test script to load exif info into Matlab
%clear all

fname = 'IR_90781.jpg';

exif = which('exiftool.exe');

TS=[ '"' exif '" -b -RawThermalImage "' fname '"']; 
[status, exifdata] = system(TS); 

% %iPad image sample
% R1 = 1931.045;
% R2 = 0.013960212;
% B = 1475.5;
% F = 1;
% O = -2535;
% S =  ;
% T_refl = 22+273.15;
% Emissivity = 0.97;

%Camera image sample
R1 = 20040.059;
R2 = 0.012000696;
B = 1485.4;
F = 1;
O = -6670;
S =  ;
T_refl = 17+273.15;
Emissivity = 0.97;
 
RAW_refl = (R1/(R2*(xp(B/T_refl)-F))- O;
RAW_obj = (S-(1-Emissivity)*RAW_refl)/Emissivity;
T_obj = B/ln((R1/R2*(RAW_obj+O))+F);
