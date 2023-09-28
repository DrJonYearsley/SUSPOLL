# SUSPOLL (Bee Model)

This directory contains a simplified rewrite of the original bee model

bee_temp.m
  The primary script that runs the model

set_parameters.m
  This is called by bee_temp.m. It sets default parameters for the model
  based upon options specified at the start of bee_temp.m
  (e.g. behaviour, species)

dTth_dt.m
  This is called by bee_temp.m via the ode solver.
  The script contains a function that calculates the rate of change
  of thorax temperature

myEvent.m
  This script is identifies when throax temperature exceeds an upper limit
  and stops the integration routine in bee_temp.m
  Using this event identifier speeds up the computations when running
  many iterattions