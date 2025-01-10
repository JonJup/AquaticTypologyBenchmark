# —————————————————————————— #
# ——— PARENT SCRIPT WP 1 ——— # 
# —————————————————————————— #

# Jonathan Jupke (jonjup@protonmail.com)
# 11.12.2024

# calls scripts ---------------------------------------------------------------------
source("code/00_setup.R", local = FALSE)

# CALL SCRIPTS ----------------------------------------------------------------------
#1. add_eu_hydro_id_to_fish
#2. combine_env 
#3. add env_fish - Adding environmental variables to fish observations 
#4. HMSC
source("code/04_HMSC.R")
#5. evaluate the likelihood of each parameter combination
source("code/05_estimate_density_w_dirichlet.R")
#5. evaluate typoloties 
source("code/06_evaluate typology.R")
#- 07. Combine Dirichlet densities with evaluation results
source("code/07_combine_dirichlet_density_and_evaluation.R")
#- PLOT 
source("code/p_01_vis_results.r")