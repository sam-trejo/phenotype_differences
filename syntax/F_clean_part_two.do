
use "${analytic_temp_plus_meta}", clear

***********************************************************************************
*** FIRST DIFFERENCES WITHIN TRANSORMATION ON OUTCOME VARIABLES
***********************************************************************************

***phenotypes
foreach var of global stub {
	***generate family N and mean of variable
	egen fam_n_`var'=count(out_`var'), by(idpub)
	egen fam_mn_`var'=mean(out_`var'), by(idpub)
	
	***generate phenotype difference variable
	gen diff_out_`var'=(out_`var'-fam_mn_`var')*2
	
	***set phenotype difference as missing if there's only one value in a family
	replace diff_out_`var'=. if fam_n_`var'!=2

	***drop extra variables
	drop fam_n_`var' fam_mn_`var'
}

*** mortality outcomes & covariates
foreach var of varlist female* c_born* imp_span alive_by_* { // pgi*
	egen fam_n_`var'=count(`var'), by(idpub)
	egen fam_mn_`var'=mean(`var'), by(idpub)

	gen diff_`var'=(`var'-fam_mn_`var')*2

	replace diff_`var'=. if fam_n_`var'!=2

	drop fam_n_`var' fam_mn_`var'
}

***********************************************************************************
*** RESIDUALIZE SIBLING OUTCOME DIFFERENCES ON COVARIATE DIFFERENCES
***********************************************************************************

***mortality outcomes
foreach var of varlist diff_imp_span diff_alive_by* {
	reg `var' diff_female diff_c_born diff_c_born2 diff_female_c_born diff_female_c_born2 if dna_sib==0
	predict res_`var' if e(sample), residuals
}

***********************************************************************************
*** RESTRICT AND SAVE OUT DATA
***********************************************************************************

***keep only relevant variables
keep individ familypub personid idpub rtype is_grad is_sib dna* srel* /// id vars
	 *died *dead age2020 span imp_span res_diff_imp_span alive* res_diff_alive* /// death vars
	 female* born c_born c_born2 /// covariates vars
	 ses57 pop57 /// milwaukee57 farm57  lnparinc57_mi parocc57_mi daded momed ocf357 pi5760
	 *out* pgi* diff* // res*  pc*
	 
***order variables
order individ familypub personid idpub rtype is_grad is_sib dna* srel* /// id vars
	  *died *dead age2020 span imp_span res_diff_imp_span alive* res_diff_alive* /// death vars
	  female* born c_born c_born2 /// covariates vars
	  ses57 pop57 /// milwaukee57 farm57  lnparinc57_mi parocc57_mi daded momed ocf357 pi5760
	  *out* *pgi* diff* //  res*  pc*  
	 
***********************************************************************************
*** MERGE ON RHO DATA
***********************************************************************************

merge 1:1 idpub is_grad using "${rho}", nogenerate keep(1 3)

***make all observations have non-missing rho for all PGIs
foreach var of varlist rho* {
	display "`var'"
	sum `var'
	replace `var' = `r(mean)'
}

***********************************************************************************
*** GENERATE RHO-TRANSFORMED PGI FOR PHENOTYPE DIFFERENCES
***********************************************************************************

foreach var of varlist pgi* {
	gen x5_`var'=(1 - rho_`var')*`var'
}

***********************************************************************************
*** SEX-SPECIFIC PGIS
***********************************************************************************

replace pgi_menses = . if female == 0
replace pgi_neb_fem = . if female == 0
replace pgi_deep = . if female == 1
replace pgi_neb_male = . if female == 1
	
***********************************************************************************
*** MERGE ON MEASUREMENT ERROR MULTIPLIER
***********************************************************************************

gen id=1
merge m:1 id using "${data}/me_mult_wide_2024_01_05.dta", nogenerate
drop id

***********************************************************************************
*** LABEL VARIABLES
***********************************************************************************
 	  
label var out_bmi "Body Mass Index"
label var out_hgt "Height"
label var out_cog "Cognitive Ability"
label var out_edu "Years of Schooling"
label var out_birth "Age at First Birth"
label var out_dep "Depressive Symptoms"
label var out_extra "Extroversion"
label var out_neuro "Neuroticism"
label var out_open "Openness to Experience"
label var out_risk "Risk Tolerance"

label var is_grad "Graduate"
label var born "Birth Year" 
label var female "Female"
label var dead "Deceased Prior to 2020"
label var alive_by_75 "Survived to Age 75"
label var imp_span "Lifespan*"
label var dna_sib "Two Genotype Sample"

***********************************************************************************
*** SAVE OUT DATA
***********************************************************************************

save "${data}/clean/analytic_${date}.dta", replace
