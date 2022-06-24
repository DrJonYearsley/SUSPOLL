
%BB thorax mass [0.014,0.132]
%HB thorax mass [0.0349,0.0465]
masses = 0.01:0.0005:0.16;   %full range
%masses = 0.001:0.05:0.2;    %testing the code limited range
pcolors=['b' 'y' 'r'];
tspan = 0:2000; 
MaxAirTemp_BB = zeros(length(masses),3);  %row for each mass, column for each flight speed
MaxAirTemp_HB = zeros(length(masses),3);  %row for each mass, column for each flight speed
BB = [1,0];  %do BB first
HB = [0,1];  %do HB second

for Bee = 1:2
%%%% Set model switches
Bumblebee = BB(Bee);    %only one of Bumblebee and Honeybee can be true
Honeybee = HB(Bee);
CoolingOff = false;   %any combination of cooling can be true
CoolingSwitch = false;  
CoolingOn = true; 
Resting = false;      %any combination of metabolic states can be true 
Shivering = false; 
Flying = true;
VaryMass = true;     %only true if varying mass of bee
Switches = [Bumblebee, Honeybee, CoolingOff, CoolingSwitch,CoolingOn,Resting,Shivering,Flying,VaryMass];

    
for m = 1:length(masses) 
    m_thorax = masses(m)  %print out to show where we are
    
    for v_bee = 1:3
    %%%% Set parameters 
    Parameter_Values = Parameters(Switches,m_thorax,v_bee);  
    
    
    %%%% Run Lookup Table
    TemperaturesTable = LookupTable_Function(Switches,Parameter_Values,tspan);
    
    
    %%%% Calculate when bee goes into thermal danger zone even with cooling
    LethalTemp = Parameter_Values(53);  %lethal thorax/air temp in C
    CoolingTemp = Parameter_Values(54);   %thorax temp where cooling begins
    FlyingTemp = Parameter_Values(55);      %thorax temp where flight can begin
    maxy=Parameter_Values(24)-273.15;    %don't solve above 60C because it's not biologically relevant

    MaxAirTemp_indices = find(TemperaturesTable(:,9) <= CoolingTemp);  %get the array indices where flying bee with cooling thorax temp >= lethal
    length_of_indices = length(MaxAirTemp_indices);
    
    if Honeybee==true 
        if isempty(MaxAirTemp_indices)
            MaxAirTemp_HB(m,v_bee) = 0;  %if it always reaches lethal point, set to min air temp
        elseif length_of_indices == 51  %all of the air temps are less than lethal thorax point
            MaxAirTemp_HB(m,v_bee) = maxy;  %if it never reaches lethal point, set to max air temp
        else
            MaxAirTemp_HB(m,v_bee) = TemperaturesTable(MaxAirTemp_indices(length_of_indices),10);  %if it does, get the air temp for the last index 
        end
    end
    
    if Bumblebee==true
        if isempty(MaxAirTemp_indices)
            MaxAirTemp_BB(m,v_bee) = 0;  %if it always reaches lethal point, set to min air temp
        elseif length_of_indices == 51  %all of the air temps are less than lethal thorax point
            MaxAirTemp_BB(m,v_bee) = maxy;  %if it never reaches lethal point, set to max air temp
        else
            MaxAirTemp_BB(m,v_bee) = TemperaturesTable(MaxAirTemp_indices(length_of_indices),10);  %if it does, get the air temp for the last index 
        end
    end
    
    end
end
end
writematrix(MaxAirTemp_BB,'MaxAirTemp_BB.csv');
writematrix(MaxAirTemp_HB,'MaxAirTemp_HB.csv');


%plot BB curves, smothed
plot(masses,smooth(MaxAirTemp_BB(:,1),0.1,'loess'),'r')
hold on
plot(masses,smooth(MaxAirTemp_BB(:,2),0.1,'loess'),'r')
plot(masses,smooth(MaxAirTemp_BB(:,3),0.1,'loess'),'r')

%plot HB curves, smothed
plot(masses,smooth(MaxAirTemp_HB(:,1),0.1,'loess'),'b')
plot(masses,smooth(MaxAirTemp_HB(:,2),0.1,'loess'),'b')
plot(masses,smooth(MaxAirTemp_HB(:,3),0.1,'loess'),'b')

%Fix window and labels
xlim([0.01 0.16])
xlabel('mass of thorax (g)')
ylabel('Maximum air temperature for sustained flying')


%add ranges to top
% y = 44;
% line([0.069*(1/2.95),0.238*(1/2.95)],[y,y], 'color', 'r')
% y = 43;
% line([0.0336,0.0485],[y,y], 'color', 'b')

x = [0.052];  %midpoint of BB range
y = [44];    %suitable y-value
err = [0.0286];  %half width of BB range
errorbar(x,y,err,'horizontal', 'color', 'r')
plot(0.057,44,'color','b','Marker','.', 'MarkerSize', 12)  %default BB m_th

x = [0.0411];  %midpoint of HB range
y = [43.5];    %suitable y-value
err = [0.0075];  %half width of HB range
errorbar(x,y,err,'horizontal', 'color', 'b')
plot(0.0407,43.5,'color','r','Marker','.', 'MarkerSize', 12)  %default HB m_th

% xline(0.014,'r-') %Bumblebee range
% xline(0.132,'r-')
% xline(0.069*(1/2.95),'r-') %Bumblebee range
% xline(0.238*(1/2.95),'r-')
%body mass 69-238mg from Heinrich1983
%scale by body mass/thorax mass = 2.95 from Kammer1974

% xline(0.0336,'b-')  %Honeybee range, thorax masses observed in Cooper1985
% xline(0.0485,'b-')