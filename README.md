Analytic method presented in "Quantifying viral pandemic potential from experimental transmission studies", currently a pre-print on bioRxiv at https://www.biorxiv.org/content/10.1101/2025.03.24.645081v1. 

In this manuscript, we aimed to extract more information from transmission experiments by (1) mathematically relating within-host viral titers in index animals to the probability of and timing of transmission to paired contact animals, (2) quanitfying "transmissibility", (3) simulating onward transmission from "naturally infected" individuals and estimating basic epidemiological parameters, and (4) extending these analyses across a range of contact rates to estimate the dynamics of an invading strain. We apply this approach to two influenza virus transmission experiments, one using A/California/07/2009 (H1N1) and one using A/California/07/2009 (H3N2). The raw data from these experiments are available as "H1N1_raw_titer_data.csv" and "H3N2_raw_titer_data.csv", respectively. 

The R files "H1N1.titers" and "H3N2.titers" visualize this data and are shown in the manuscript as Figure 1. 

The R file "functional_forms" tests the linear, threshold, and Hill functional forms, and visualizes their relationship to the probability of transmission for given parameterizations, shown in Figure 3.

The R file "MLE" runs a maximum likelihood estimation for "s" (the log functional form) for both H1N1 and H3N2. This is shown in the manuscript as Figure 2.

The R file "assessing_pandemic_potential" runs the analyses shown in Figure 4. 

The "Epi_characteristics_experiment" folder contains the R files "varying_contact_sims" and "invasion_dynamics". The "varying_contact_sims" runs transmission simulations (equivalent to those in "assessing_pandemic_potential") at multiple contact rates. The "invasion_dynamics" file runs the analyses shown in Figure 5, and uses the results from the varying contact simulations.

The "Supplemental" folder contains R files involved in analyses shown in the supplement. In particular, "contact_dynamics" shows the results from Figure S1, while "func_forms_pandemic_potential" and "func_forms_contact_sims" run transmission simulations equivalent to "assessing_pandemic_potential" and "varying_contact_sims" for the linear, threshold, and Hill functional forms.
