%DesiredTemps = [10 25 35];  %temperatures I'd like to match in model 
%air temps at which [shivering bee can warm to 30C, abdomen cooling starts, bee goes above 42C]

% %note: no data for third air temp point, so leave it out for now
% DesiredTemps = [6 25];  %temperatures I'd like to match in model 

s_default = 0.9;  %from Church1960
I_flying = 0.06229515;     %Kammer1974, converted to W

%grid for i0 and f
% i0_min = 0.001349728;  %shouldn't go lower than resting i0... 
% i0_max = 0.15;  %Heinrich number
i0_min = 0.001;  
i0_max = 0.01;  

s_min = 0.85;
s_max = 1;

i0_axis = [i0_min:0.00005:i0_max];
s_axis = [s_min:0.0005:s_max];

i0_grid_size = length(i0_axis);
s_grid_size = length(s_axis);

%Storage matrices
%rows for i0, columns for r
ShiverTo30 = zeros(i0_grid_size,s_grid_size);  
Diverge = zeros(i0_grid_size,s_grid_size);
ThermalDanger = zeros(i0_grid_size,s_grid_size);

%tic
%ticBytes(gcp);
parfor i = 1:i0_grid_size
    for s = 1:s_grid_size
        %[i,f]
        i0_guess = i0_axis(i);
        s_guess = s_axis(s);
        KeyTemps = RunModelABC(i0_guess,s_guess);
        ShiverTo30(i,s) = KeyTemps(1); %how far from warming goal
        Diverge(i,s) = KeyTemps(2);  %how far from divergence goal
        ThermalDanger(i,s) = KeyTemps(3);  %temp where bee goes above 42C
    end
end
%tocBytes(gcp);
%toc
writematrix(ShiverTo30,'ShiverTo30_i0_s_2.csv');
writematrix(Diverge,'Diverge_i0_s_2.csv');
writematrix(ThermalDanger,'ThermalDanger_i0_s_2.csv');

% ShiverTo30_dist = abs(ShiverTo30-DesiredTemps(1));
% Diverge_dist = abs(Diverge-DesiredTemps(2));
% writematrix(ShiverTo30_dist,'ShiverTo30_dist_i0_s_1.csv');
% writematrix(Diverge_dist,'Diverge_dist_i0_s_1.csv');
