# Meritocracy
Data and code for "Theoretical foundation and empirical assessment of surname representation and meritocracy in academia"
EW Tekwa, Rachel Giles, Alexandra Davis
Sep 7, 2022

- Raw data was collected with source links and date in "Surname citation_full.xlsx"
- To run simulation models, run the Matlab script "runNameRepresentation_factorial.m". This script calls "nameRepresentation_noPlot.m" or "nameRepresentation_noPlot_lognorm.m" which runs each simulation replicate. Change parameter values in lines 12-54 for sensitivity tests.
- All data from the paper's main result simulations is contained in "simulations_mainResults.mat". This file can be loaded into the Matlab workspace without running new simulations to explore existing results. Run line 200 onward in "runNameRepresentation_factorial.m" to plot results.
- To run empirical analyses, run the Matlab script "runEmpiricalSurnameStatsStack.m". Change line 38 to plot LR data with different reference population normalizations. Script automatically loads "EmpiricalSurnameData.mat" that contains all empirical data analyzed in the paper.
