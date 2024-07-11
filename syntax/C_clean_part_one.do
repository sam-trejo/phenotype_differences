/*************************************************************************
Goal: This code creates a clean WLS data file for phenotype differences analysis. 
**************************************************************************/

***********************************************************************************
*** LOAD WLS DATA
***********************************************************************************

***load main wls data
use "${data}\wls\wls_bl_14.03.stata\wls_bl_14_03.dta", clear

***********************************************************************************
*** MERGE ON AUXILLARY DATA
***********************************************************************************

capture drop _merge

***merge on pgi repo
merge m:1 idpub rtype using "${data}\pgs\PGIrepo_idpub_v1.1\PGIrepo_v1.1_idpub_shuffled.dta", keep(1 3) nogenerate

***merge on Klint's phenotype data
merge 1:1 idpub rtype personid using "${data}\clean\wls_phenotypes.dta", keep(3) nogenerate

***********************************************************************************
*** CLEAN VARIABLES
***********************************************************************************

egen individ = group(idpub rtype)
isid individ

***grad vs. sib
gen is_grad = rtype=="g"
gen is_sib = rtype=="s"

***birth year
generate born = z_brdxdy
drop if born==.

***centered birth year (to distribution of graduates)
gen c_born = .
sum born if rtype=="g"
replace c_born = (born - `r(mean)')

***sex
generate female = .
replace female = 0 if z_sexrsp==1
replace female = 1 if z_sexrsp==2
drop if female==.

***sex and birth year interactions
generate c_born2 = c_born * c_born
generate female_c_born = female * c_born 
generate female_c_born2 = female * c_born2

***death date
gen dead = z_livgrad - 1
gen died = z_deatyr if inrange(z_deatyr,1900,2100)
replace died = 2019 if inrange(died,2020, 2022)

***lifespan
gen span = died - born

***age in 2020 (even for deceased folks)
gen age2020 = 2020 - born

***********************************************************************************
*** 
***********************************************************************************

***standardize pgi and outcomes (to distribution of graduates)
foreach var of varlist pgi* out* {
	sum `var' if rtype=="g"
	replace `var' = (`var' - `r(mean)')/`r(sd)'
}

***generate dna sample indicator
gen dna = 1
foreach var of varlist pgi* {
    replace dna = 0 if `var'==.
}

***save out full sample for selection into genotyping
preserve
	drop  z_brdxdy - z_q1a203re
	save "${data}\clean\full_${date}.dta", replace
restore

***keep only first two people in each family
tab personid
keep if personid<=2

***keep only self-reported full sibling pairs
keep if srel1==0 | srel1==1
keep if srel2==0 | srel2==1

egen count=count(individ), by(idpub)
keep if count==2
drop count

******create indicator for sibling pairs that are both genotyped
egen total=total(dna), by(idpub)
tab total
keep if inlist(total,1,2)
gen dna_sib = 0
replace dna_sib = 1 if total==2
drop total

***create indicator for sibling pairs that have experiencd at least one death
egen sib_tot_dead = total(dead), by(idpub)
gen sib_died = 0
replace sib_died = 1 if inlist(sib_tot_dead,1,2)
drop sib_tot_dead

***********************************************************************************
*** 
***********************************************************************************

***either lifespan or age in 2022
gen imp_span = span
replace imp_span = age2020 if died==.



***generate alive by 30-80
forval i=30(5)80 {
	gen alive_by_`i' = span>=`i'
	replace alive_by_`i' = . if age2020<`i' & (dead==0)
}

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

***********************************************************************************
*** 
***********************************************************************************

save "${data}\clean\analytic_temp_${date}.dta", replace
