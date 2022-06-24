clear all
%T= B/ln(R1/(R2*(S+O))+F)
% T= object temperature in Kelvins
% S = 16 Bit RAW value
% R1 = Planck R1 constant
% R2 = Planck R2 constant
% B = Planck B constant. Value range 1300-1600
% F = Plank F constant. Value range 0.5 - 2.
% O = Planck O (offset) constant. Its a negative value
% ln() = natural logarithm

%system('/usr/bin/exiftool -PlanckR2 IR_7701.jpg')

imagename = 'IR_90812';


exifpath = 'C:\Users\UCD\Documents\SUSPOLL\Fieldwork\PilotThermalCamera\Camera\exiftool.exe';
dir = 'C:\Users\UCD\Documents\SUSPOLL\Fieldwork\PilotThermalCamera\Camera\2021-06-19';
% exifpath = 'C:\Users\UCD\Documents\SUSPOLL\Fieldwork\PilotThermalCamera\iPad\exiftool.exe';
% dir = 'C:\Users\UCD\Documents\SUSPOLL\Fieldwork\PilotThermalCamera\iPad';
sep='\';

fname = [dir  sep imagename '.jpg'];
fname2 = [ dir sep imagename '.tif'];

str = ['exiftool -b -RawThermalImage ' dir sep imagename '.jpg > ' dir sep imagename '.tif'];
system(str);

% Read in parameters from the jpg file
param = {'R1','R2','B','F','O'};
for p=1:length(param) 
    str1 = ['"' exifpath '" -Planck' param{p} ' "' fname '" '];
    [a,b ] = system(str1);
    tmp = strsplit(b,':');
    eval([param{p} '= str2double(tmp{2});']);
end

str2 = ['"' exifpath '" -ReflectedApparentTemperature' ' "' fname '" '];
[c,d ] = system(str2);
tmp2 = strsplit(d);
T_refl = str2double(tmp2{5});
%% manually get all the values and calculate temperature
Em = 0.97;  %real emmissivity of bees, not the camera default
S_tmp=imread(fname2);
S = S_tmp(213,317);
RAW_refl = R1/(R2*(exp(B./T_refl)-F)) - O;
RAW_obj = (double(S_tmp) - (1-Em)*RAW_refl)./Em;
T = real(B ./ log(R1./(R2 *(RAW_obj + O))+F)-273.15);
%T = real(B ./ log(R1./(R2 *(RAW_obj + O))+F));


imagesc(T)
colorbar()
% caxis([25,45])

% impixelinfo



% impixelinfo
X = T(240:285,290:330);  %get the thorax - have to identify by hand
imagesc(X)
[m,n] = size(X); 
TN = m*n;  %number of pixels in it
TTP = TN*.1;  %10% of those pixels
TTP = ceil(TTP); %round in case TN is an odd number and TTP is not an integer

VX = X(:);               %places values in single column
[val,ind] = sort(VX, 'descend');   %sort them in decreasing order
y = VX(ind(1:TTP));  %take the top 10%
mean(y)  %average temp
writematrix(VX,[ dir sep imagename '.csv']);





%T=real(1514./log(19242.684./(0.076735988*(double(S)-7722))+1.5))-272.15;
% convert Kelvin to celsius
%T2=T-272.15;

%% read in jpg
% S2 = imread(fname);
% 
% T2 = dlmread('/Users/jon/Dropbox/smouldering/data/burn10_IR_7701.csv',',',6,0);
% 
% %% add temperatures from .jpg image
% % import .csv file = variable b
% 
% 
% h = 6.62606957E-34;
% kb = 1.3806488E-23;
% lambda = 9.5E-6;
% c = 299792458;
% B2 = h*c/(lambda*kb);
% 
% 2*h*c^2/lambda^5;
