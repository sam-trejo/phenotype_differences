***********************************************************************************
*** IMPORT KLINT PHENOTYPE DATA
***********************************************************************************

import delimited "${data}\clean\wls_phenotypes.csv", clear stringcols(_all) 

gen neb_male = number_ever_born if female=="0"

gen neb_fem = number_ever_born if female=="1"

drop v1 id female number_ever_born

***********************************************************************************
*** RENAME AND LABEL VARS
***********************************************************************************

rename cognitive_performance cog
rename ed_attainment edu
rename age_first_birth birth
rename age_first_menses menses
rename alcohol_misuse alcoh
rename allergy_cat cat
rename allergy_dust dust
rename allergy_pollen pollen
rename hayfever hay
rename height hgt
rename asthma asthma
rename asthma_eczema_rhinitis aer
rename cigarettes_per_day cig_day
rename depressive_symptoms_alt dep
rename drinks_per_week drinks
rename ever_smoker ever_smk
rename physical_activity phys_act
rename self_rated_health health
rename extraversion extra
rename neuroticism neuro
rename openness open
rename life_satisfaction_family sat_fam
rename life_satisfaction_finance sat_fin
rename life_satisfaction_job sat_job
rename loneliness lonely
rename religious_attendance relig
rename risk_tolerance risk
rename subjective_well_being swb

foreach var of varlist bmi-neb_fem {
	label var `var' "`var'"
	rename `var' out_`var'
}

rename out_depressive_symptoms depressive_symptoms

label var out_cig_day "cig-day"
label var out_ever_smk "ever-smk"
label var out_phys_act "phys-act"
label var out_neb_male "neb-male"
label var out_neb_fem "neb-fem"
label var out_sat_fam "sat-fam"
label var out_sat_fin "sat-fin"
label var out_sat_job "sat-job"

***********************************************************************************
*** DESTRING AND SAVE OUT
***********************************************************************************

destring idpub familypub personid out_bmi-out_neb_fem, force replace

replace out_sat_job = -1 * out_sat_job

save "${data}\clean\wls_phenotypes.dta", replace

