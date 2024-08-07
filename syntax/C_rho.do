***CORRELATION TABLE

***********************************************************************************
*** LOAD WLS AND PGI DATA 
***********************************************************************************

matrix drop _all

***load main wls data
use "${data}/wls/wls_bl_14.03.stata/wls_bl_14_03.dta", clear

***merge on pgi repo
merge m:1 idpub rtype using "${data}/pgs/PGIrepo_idpub_v1.1/PGIrepo_v1.1_idpub_shuffled.dta", keep(3) nogenerate

***generate WLS graduate indicator
gen is_grad = rtype=="g"

***restrict variables
keep idpub is_grad pgi*
order idpub is_grad pgi*

***********************************************************************************
*** RESTRICT TO ONE PGI PER TRAIT & RENAME TO MATCH OUTCOMES
***********************************************************************************

drop pgi_adhdsingle pgi_adventuresingle pgi_afbsingle pgi_asteczrhisingle pgi_asthmasingle ///
	 pgi_auditsingle pgi_cpsingle pgi_depsingle pgi_dpwsingle pgi_extrasingle /// 
	 pgi_famsatsingle pgi_friendsatsingle pgi_hayfeversingle pgi_highmathsingle pgi_leftoutsingle ///
	 pgi_menarchesingle pgi_nebwomensingle pgi_neurosingle pgi_religattsingle pgi_risksingle ///
	 pgi_selfhealthsingle pgi_selfmathsingle pgi_swbsingle pgi_eamulti

/*
keep pgi_activitysingle pgi_bmisingle pgi_cannabissingle pgi_cpdsingle pgi_eversmokesingle /// pgi vars
	 pgi_heightsingle pgi_migrainesingle pgi_morningsingle pgi_narcissingle pgi_nearsightedsingle /// pgi vars
 	 pgi_opensingle pgi_readingsingle pgi_easingle pgi_adhdmulti pgi_adventuremulti
	 pgi_afbmulti pgi_allergycatmulti pgi_allergydustmulti pgi_allergypollenmulti pgi_asteczrhimulti ///
	 pgi_asthmamulti pgi_auditmulti pgi_cogempmulti pgi_copdmulti pgi_cpmulti ///
	 pgi_delaydiscmulti pgi_depmulti pgi_dpwmulti pgi_extramulti pgi_famsatmulti /// 
	 pgi_finsatmulti pgi_friendsatmulti pgi_hayfevermulti pgi_highmathmulti pgi_leftoutmulti ///
	 pgi_lonelymulti pgi_menarchemulti pgi_nebmenmulti pgi_nebwomenmulti pgi_neuromulti ///
	 pgi_religattmulti pgi_riskmulti pgi_selfhealthmulti pgi_selfmathmulti pgi_swbmulti ///
	 pgi_voicedeepmulti pgi_worksatmulti
*/	  

***rename PGIs w/phenotype data
rename pgi_bmisingle pgi_bmi
rename pgi_cpmulti pgi_cog 
rename pgi_easingle pgi_edu 
rename pgi_afbmulti pgi_birth 
rename pgi_menarchemulti pgi_menses 
rename pgi_nebmenmulti pgi_neb_male
rename pgi_nebwomenmulti pgi_neb_fem
rename pgi_auditmulti pgi_alcoh
rename pgi_allergycatmulti pgi_cat 
rename pgi_allergydustmulti pgi_dust 
rename pgi_allergypollenmulti pgi_pollen
rename pgi_hayfevermulti pgi_hay
rename pgi_asthmamulti pgi_asthma
rename pgi_asteczrhimulti pgi_aer
rename pgi_cpdsingle pgi_cig_day
rename pgi_copdmulti pgi_copd
rename pgi_dpwmulti pgi_drinks
rename pgi_eversmokesingle pgi_ever_smk
rename pgi_activitysingle pgi_phys_act
rename pgi_selfhealthmulti pgi_health
rename pgi_extramulti pgi_extra
rename pgi_neuromulti pgi_neuro
rename pgi_opensingle pgi_open
rename pgi_finsatmulti pgi_sat_fam
rename pgi_famsatmulti pgi_sat_fin
rename pgi_worksatmulti pgi_sat_job
rename pgi_lonelymulti pgi_lonely
rename pgi_religattmulti pgi_relig 
rename pgi_riskmulti pgi_risk
rename pgi_swbmulti pgi_swb 
rename pgi_heightsingle pgi_hgt
rename pgi_depmulti pgi_dep

***rename PGIs w/o phenotype data
rename pgi_adhdmulti pgi_adhd
rename pgi_cannabissingle pgi_canna
rename pgi_migrainesingle pgi_migrn
rename pgi_morningsingle pgi_chrono
rename pgi_narcissingle pgi_narci
rename pgi_nearsightedsingle pgi_near_sgt
rename pgi_readingsingle pgi_read
rename pgi_adventuremulti pgi_adv
rename pgi_cogempmulti pgi_cog_emp
rename pgi_delaydiscmulti pgi_dly_disc
rename pgi_friendsatmulti pgi_sat_frnd
rename pgi_highmathmulti pgi_high_math
rename pgi_leftoutmulti pgi_leftout
rename pgi_selfmathmulti pgi_self_math
rename pgi_voicedeepmulti pgi_deep	  

***hold data
tempfile hold
save `hold', replace

***********************************************************************************
*** CREATE LIST OF PHENOTYPES STRING
***********************************************************************************

global pgi_name ""
global pgi_num ""
global comma ""

local i=1

foreach var of varlist pgi* {
	
	label var `var' "`var'"
	display "`var'"
	global pgi_name "$pgi_name `var'"
 	global comma "$comma `var',"
	
	***add leading 0 to 1-9 for variable names
	local ilead = string(`i',"%02.0f")
	
	***rename polygenic scores
	rename `var' pgi`ilead'_
	
	global pgi_num "$pgi_num pgi`ilead'_"
	
	local i = `i' + 1
}

display "$pgi_name"
display "$pgi_num"
display "$comma"

global comma = substr("$comma", 1, length("${comma}")-2)

display "$comma"

***********************************************************************************
*** RESHAPE WIDE SO EACH ROW IS A SIBLING PAIR
***********************************************************************************

***reshape data so each row is a sibling pair
reshape wide pgi*, i(idpub) j(is_grad)

***********************************************************************************
*** RANDOMLY ASSIGN SIBLING ONE AND SIBLING TWO AND ESTIMATE CORRELATION 1000 TIMES
***********************************************************************************

***
gen y=.
gen x=.
foreach var in $pgi_num {
	
	***
	center `var'0 `var'1, standardize inplace
	
	***
	forval i=1/${corr_reps} { 
		quietly {
			***
			gen rand=rnormal()
			egen mn_rand=mean(rand), by(idpub)
			replace x=`var'0 if mn_rand<0
			replace y=`var'1 if mn_rand<0
			replace x=`var'1 if mn_rand>=0
			replace y=`var'0 if mn_rand>=0
			drop *rand
			
			***
			reg y x
			local b=_b[x]
			local se=(_se[x])^2
			matrix `var'=nullmat(`var') \ `b', `se' 
		}
	}
}

***********************************************************************************
*** COMBINES ESTIMATES
***********************************************************************************

***combine all phenotypes into one big matrix
foreach var in $pgi_num {
	mat combine = nullmat(combine), `var' 
}

***average sibling correlations within phenotype
mat row = J(rowsof(combine),1,1)
mat col = row'*combine
mat corr=col/rowsof(combine)

***convert variance to standard error
forvalues i=1/47 {
	mat out = nullmat(out), corr[1,2*`i'-1], corr[1,2*`i']^.5
}

***********************************************************************************
*** EXPORT CORRELATION DATA
***********************************************************************************

use `hold', clear

***convert rhos from matrix to variables
svmat corr

foreach var of varlist corr* {
	sum `var'
	replace `var' = `r(mean)' if `var'==.
}

local i=0
foreach var of varlist pgi* {
	
	***reattach phenotype names 
	local i = `i' + 1
	rename corr`i' rho_`var'
	
	local i = `i' + 1
	rename corr`i' var_rho_`var'
	gen se_rho_`var' = var_rho_`var'^.5
}

drop var*

save "${data}/clean/phen_diff_rho_${date}.dta", replace
