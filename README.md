# Meritocracy
Data and code for "Theoretical foundation and empirical assessment of surname representation and meritocracy in academia"
EW Tekwa, Rachel Giles, Alexandra Davis
Sep 7, 2022

- Raw data was collected with source links and date in "Surname citation_full.xlsx"
- To run simulation models, run the Matlab script "runNameRepresentation_factorial.m". This script calls "nameRepresentation_noPlot.m" or "nameRepresentation_noPlot_lognorm.m" which runs each simulation replicate. Change parameter values in lines 12-54 for sensitivity tests.
- All data from the paper's main result simulations is contained in "simulations_mainResults.mat". This file can be loaded into the Matlab workspace without running new simulations to explore existing results. Run line 200 onward in "runNameRepresentation_factorial.m" to plot results.
- All sensitivity test simulations are contained in the following files that can be loaded and plotted as in the previous step:
  - smaller population: "1944 simulations_allStats_t20_numNames200_AcadPorp02_meritCapCausal_full.mat"
  - surname inheritance by capital: "1944 simulations_allStats_t20_numNames1000_AcadPorp02_meritCapCausal_CapSurname.mat"
  - lognormal capital and merit distributions: "1944 simulations_allStats_t20_numNames1000_AcadPorp02_meritCapCausal_lognorm.mat"
  - measurements at 5th generation: "1944 simulations_allStats_t5_numNames1000_AcadPorp02_meritCapCausal.mat"
  - measurements at 15th generation: "1944 simulations_allStats_t15_numNames1000_AcadPorp02_meritCapCausal.mat"
  - measurements at 25th generation: "1944 simulations_allStats_t25_numNames1000_AcadPorp02_meritCapCausal.mat"
- To run empirical analyses, run the Matlab script "runEmpiricalSurnameStatsStack.m". Change line 38 to plot LR data with different reference population normalizations. Script automatically loads "EmpiricalSurnameData.mat" that contains all empirical data analyzed in the paper.
