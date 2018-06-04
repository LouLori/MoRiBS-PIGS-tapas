# MoRiBS-PIGS
1. I have added three preprocessors: VPOTTWOLINEARROTORS, VPOTTWOLINEARROTORSIO and GETPOT. 

2. Here I have corrected the unit of potential energy (cm-1 to Kelvin) by dividing "potl" in mc_estim.cc and mc_piqmc.cc by "Units.kelvin". I have defined new unit transformer "Units.kelvin" for cm-1 to Kelvin.  In vh2h2.f, unit of potl is cm-1. 


3. GETPOT preprocessor give us to generate potential and total energy as a function "Inter molecular distances". The "Inter molecular distances" is included by qmc.input. 

4. For job submission and analyzing output data (pigs.eng), I have written two script files in python: script-for-analysis.py and script-for-job-submission.py.
