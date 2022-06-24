%DesiredTemps = [10 25 35];  %temperatures I'd like to match in model 
%air temps at which [shivering bee can warm to 30C, abdomen cooling starts, bee goes above 42C]

%note: no data for third air temp point, so leave it out for now
DesiredTemps = [6 25];  %temperatures I'd like to match in model 

r_default = 0.0367/9;  %0.0041 Calculated from Heinrich1976, 2.2 J/min at T_th-T_abdomen = 9C
I_flying = 0.06229515;     %Kammer1974, converted to W

%grid for i0 and r
i0_min = 0.001349728;  %shouldn't go lower than resting i0... 
i0_max = 0.06229515;  %Kammer number

r_min = 0;
r_max = 0.005;

i0_axis = [i0_min:0.00001:0.0045];
r_axis = [r_min:0.000005:0.001];

i0_grid_size = length(i0_axis);
r_grid_size = length(r_axis);

%Storage matrices
%rows for i0, columns for r
ShiverTo30 = zeros(i0_grid_size,r_grid_size);  
Diverge = zeros(i0_grid_size,r_grid_size);
ThermalDanger = zeros(i0_grid_size,r_grid_size);

%tic
%ticBytes(gcp);
parfor i = 1:i0_grid_size
    for r = 1:r_grid_size
        %[i,r]
        i0_guess = i0_axis(i);
        r_guess = r_axis(r);
        KeyTemps = RunModelABC(i0_guess,r_guess);
        ShiverTo30(i,r) = KeyTemps(1); %how far from warming goal
        Diverge(i,r) = KeyTemps(2);  %how far from divergence goal
        ThermalDanger(i,r) = KeyTemps(3);  %temp where bee goes above 42C
    end
end
%tocBytes(gcp);
%toc
writematrix(ShiverTo30,'ShiverTo30_fine2.csv');
writematrix(Diverge,'Diverge_fine2.csv');
writematrix(ThermalDanger,'ThermalDanger_fine2.csv');

ShiverTo30_dist = abs(ShiverTo30-DesiredTemps(1));
Diverge_dist = abs(Diverge-DesiredTemps(2));
writematrix(ShiverTo30_dist,'ShiverTo30_dist_fine2.csv');
writematrix(Diverge_dist,'Diverge_dist_fine2.csv');

%%%%%%%%%%%%%%%%%%%%% Finer grid   %%%%%%%%%%%%%%%%%%%%%%%%%

i0_axis = [i0_min:0.0002:i0_max];
r_axis = [r_min:0.00002:r_max];

i0_grid_size = length(i0_axis);
r_grid_size = length(r_axis);

%Storage matrices
%rows for i0, columns for r
ShiverTo30 = zeros(i0_grid_size,r_grid_size);  
Diverge = zeros(i0_grid_size,r_grid_size);
ThermalDanger = zeros(i0_grid_size,r_grid_size);

%tic
%ticBytes(gcp);
parfor i = 1:i0_grid_size
    for r = 1:r_grid_size
        %[i,r]
        i0_guess = i0_axis(i);
        r_guess = r_axis(r);
        KeyTemps = RunModelABC(i0_guess,r_guess);
        ShiverTo30(i,r) = KeyTemps(1); %how far from warming goal
        Diverge(i,r) = KeyTemps(2);  %how far from divergence goal
        ThermalDanger(i,r) = KeyTemps(3);  %temp where bee goes above 42C
    end
end
%tocBytes(gcp);
%toc
writematrix(ShiverTo30,'ShiverTo30.csv');
writematrix(Diverge,'Diverge.csv');
writematrix(ThermalDanger,'ThermalDanger.csv');

ShiverTo30_dist = abs(ShiverTo30-DesiredTemps(1));
Diverge_dist = abs(Diverge-DesiredTemps(2));
writematrix(ShiverTo30_dist,'ShiverTo30_dist.csv');
writematrix(Diverge_dist,'Diverge_dist.csv');

