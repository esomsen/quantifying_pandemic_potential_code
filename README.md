Analytic method presented in "Quantifying viral pandemic potential from experimental transmission studies", currently a pre-print on bioRxiv at https://www.biorxiv.org/content/10.1101/2025.03.24.645081v1. 

In this manuscript, we aimed to extract more information from transmission experiments by (1) mathematically relating within-host viral titers in index animals to the probability of and timing of transmission to paired contact animals, (2) quanitfying "infectiousness", (3) simulating onward transmission from "naturally infected" individuals and estimating basic epidemiological parameters, and (4) extending these analyses across a range of contact rates to estimate the dynamics of an invading strain. We apply this approach to two influenza virus transmission experiments, one using A/California/07/2009 (H1N1) and one using A/California/07/2009 (H3N2). The raw data from these experiments are available as "H1N1_raw_titer_data.csv" and "H3N2_raw_titer_data.csv", respectively. 

The R files "H1N1.titers" and "H3N2.titers" visualize this data and are shown in the manuscript as Figure 1. 

The R file "MLE" runs a maximum likelihood estimation for "s" for both H1N1 and H3N2. This is shown in the manuscript as Figure 2.

The R file "assessing_pandemic_potential" runs the analyses shown in Figure 3. 

The "Epi_characteristics_experiment" folder contains the R files "generation_time_simulation" and "invasion_dynamics". The "generation_time_simulation" file estimates the generation times for H1N1 and H3N2 assuming a high contact rate, as described in the methods. The "invasion_dynamics" file runs the analyses shown in Figure 4, and uses the results from the high contact simulations for panel D (the values are hard coded). 
