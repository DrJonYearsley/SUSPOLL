
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Initialize
%DesiredTemps = [10 25 35];  %temperatures I'd like to match in model 
%air temps at which [shivering bee can warm to 30C, abdomen cooling starts, bee goes above 42C]

%note: no data for third air temp point, so leave it out for now
DesiredTemps = [6 25];  %temperatures I'd like to match in model 

T = 4;   %number of iterations after 1st one
tolerances = [8, 6, 5, 2];    %tolerence for each iteration
N = 10;  %number of samples to take in each iteration 

%prior distribution for i0
i0_min = 0;
i0_max = 0.06229515;  %Kammer number
%prior distribution for r
r_min = 0;
r_max = 0.5;

  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Run fitting... 

%for the first iteration, use uniform distributions for i0 and r
CurrentSample = zeros(N,5); 
initial_tolerance = 10;  %the initial tolerance

for i = 1:N
howmany = 0;  %to tell how many tries were needed to meet tolerance for this sample
dist_metric1 = initial_tolerance+1;  %just needs to be something greater than the initial tolerance
dist_metric2 = initial_tolerance+1;  %just needs to be something greater than the initial tolerance
    while any([dist_metric1>=initial_tolerance dist_metric2>=initial_tolerance]) %do this until the distance metric is within the tolerance
       %sample i0 and r
       i0_guess = random('Uniform',i0_min,i0_max); 
       r_guess = random('Uniform',r_min,r_max); 
       %run the model 
       KeyTemps = RunModelABC(i0_guess,r_guess);
       %dist_metric = sum((DesiredTemps - KeyTemps(1:2)).^2);
       dist_metric1 = abs(DesiredTemps(1)-KeyTemps(1));
       dist_metric2 = abs(DesiredTemps(2)-KeyTemps(2));
       howmany = howmany+1;
       [i howmany dist_metric1 dist_metric2  i0_guess r_guess] %i for sample number, howmany for how many tries
    end %end of while loop
CurrentSample(i,:) = [i0_guess, r_guess, KeyTemps]; 
end 

Sample1 = CurrentSample;

%for the rest of the iterations, sample from the previous set of samples
for t = 1:T  %iterations    
    PreviousSample = CurrentSample;   %make the old current sample the previous one
    CurrentSample = zeros(N,5);       %new current sample is empty
    dist_metric1 = tolerances(t)+2;  %just needs to be something greater than the initial tolerance
    dist_metric2 = tolerances(t)+2;  %just needs to be something greater than the initial tolerance

    for i = 1:N  %samples
        howmany = 0;  %to tell how many tries were needed to meet tolerance for this sample
        while any([dist_metric1>=tolerances(t) dist_metric2>=tolerances(t)]) %do this until the distance metric is within the tolerance
        %sample i0 and r from PreviousSample set and perturb
        i0_guess = datasample(PreviousSample(:,1),1) + random('Uniform',-0.05,0.05); 
        r_guess = datasample(PreviousSample(:,2),1) + random('Uniform',-0.05,0.05); 
        %run the model 
        KeyTemps = RunModelABC(i0_guess,r_guess);
        %dist_metric = sum((DesiredTemps - KeyTemps(1:2)).^2);
        dist_metric1 = abs(DesiredTemps(1)-KeyTemps(1));
        dist_metric2 = abs(DesiredTemps(2)-KeyTemps(2));
        howmany = howmany+1;
        [t i howmany dist_metric1 dist_metric2 i0_guess r_guess]
        end %end of while loop
        %save the sample (the one that met the tolerence)    
        CurrentSample(i,:) = [i0_guess, r_guess, KeyTemps]; 
    end %end of samples for loop
end %end of iterations for loop
