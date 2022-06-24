function cooling = Cooling_Flux(Bumblebee,r,M_th,c,A,B,T_aK,rh,Pr,MM_air,MM_vapor,R_specific,R_0,D_A,h_fg,v)

    if Bumblebee == true
        Ab = r*(1/(M_th*c)); 
        Ev1 = 0;
        Ev2 = 0;
    else  %it must be a honeybee
        Ab = 0;   %value of abdomen cooling will be 0 if Bumblebee=false

        p_sat_air = A*exp(B/T_aK);  %partial pressure of saturated vapor (humid air at saturation)
        p_air = rh*p_sat_air;   %partial pressure of humid air at temperature T_aK
        rho_humid = ((Pr-p_air)*MM_air + p_air*MM_vapor)/(R_specific*T_aK);   %density of humid air

        X_air = p_air/Pr;  %mole fraction of ambient air
        Y_air = 1/( 1 + ((1-X_air)/X_air)*(MM_air/MM_vapor) );

        m_evap_coef = 2*pi*R_0*rho_humid*D_A;
        Ev1 = h_fg*m_evap_coef*log(1-Y_air)/(M_th*c);   %constant for a fixed T_air
        Ev2 = -h_fg*m_evap_coef/(M_th*c);   %*log(1-Y_sfc)
    end     
    cooling=[Ab, Ev1, Ev2];