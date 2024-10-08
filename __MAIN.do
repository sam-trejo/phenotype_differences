/* REQUIRED PACKAGES
matselrc
frmttable
rscript
reghdfe
ftools
estout
center
*/

***********************************************************************************
*** SETUP
***********************************************************************************

clear all
set more off, perm
set maxvar 30000
capture set trace off
capture postutil close
set scheme stcolor
graph set window fontface "Calibri Light"
set seed 19146

***********************************************************************************
*** SET DIRECTORY GLOBALS
***********************************************************************************

***set files path gobals
global phen_diff "`c(pwd)'"
global data "${phen_diff}/data"
global syntax "${phen_diff}/syntax"
global figure "${phen_diff}/figures"
global table "${phen_diff}/tables"

***********************************************************************************
*** DATE GLOBAL
***********************************************************************************

*makes stata date in numeric year_month_day format
quietly {
	global date=c(current_date)

	***day
	if substr("$date",1,1)==" " {
		local val=substr("$date",2,1)
		global day=string(`val',"%02.0f")
	}
	else {
		global day=substr("$date",1,2)
	}

	***month
	if substr("$date",4,3)=="Jan" {
		global month="01"
	}
	if substr("$date",4,3)=="Feb" {
		global month="02"
	}
	if substr("$date",4,3)=="Mar" {
		global month="03"
	}
	if substr("$date",4,3)=="Apr" {
		global month="04"
	}
	if substr("$date",4,3)=="May" {
		global month="05"
	}
	if substr("$date",4,3)=="Jun" {
		global month="06"
	}
	if substr("$date",4,3)=="Jul" {
		global month="07"
	}
	if substr("$date",4,3)=="Aug" {
		global month="08"
	}
	if substr("$date",4,3)=="Sep" {
		global month="09"
	}
	if substr("$date",4,3)=="Oct" {
		global month="10"
	}
	if substr("$date",4,3)=="Nov" {
		global month="11"
	}
	if substr("$date",4,3)=="Dec" {
		global month="12"
	}

	***year
	global year=substr("$date",8,4)

	global date="$year"+"_"+"$month"+"_"+"$day"
}

dis "$date"

***********************************************************************************
*** SET DATA GLOBALS
***********************************************************************************

***set dataset globals

***
global rho "${data}/clean/phen_diff_rho_${date}.dta"

***
global analytic_temp "${data}/clean/analytic_temp_${date}.dta"
global analytic_temp_plus_meta "${data}/clean/analytic_temp_plus_meta_${date}.dta"
global analytic "${data}/clean/analytic_${date}.dta"

***
global full "${data}/clean/full_${date}.dta"

***********************************************************************************
*** SET REPETITIONS GLOBALS
***********************************************************************************

***
global corr_reps = 1000

***
global fig1C_reps = 1000

***********************************************************************************
*** SET PGI PHENOTYPES NAMES GLOBALS
***********************************************************************************

*** pgi phenotypes we observe in WLS
global stub ///
	   aer alcoh asthma birth bmi ///
	   cig_day cog copd dep drinks ///
	   edu ever_smk extra health hgt ///
	   lonely neuro open phys_act relig ///
	   sat_fam sat_fin sat_job swb
	   /// neb_male neb_fem menses [SEX SPECIFIC]
	   /// pollen risk [NO SIGNIFICANT BETWEEN-FAMILY PREDICTION]
	    // hay dust cat [2X GENOTYPE SAMPLE HAS N<300]
	 	
*** all pgi phenotypes in repository
global stub2 ///
	   adhd adv aer alcoh asthma ///
	   birth bmi canna cat chrono ///
	   cig_day cog cog_emp copd deep ///
	   dep dly_disc drinks dust edu ///
	   ever_smk extra hay health hgt ///
	   high_math leftout lonely menses migrn ///
	   narci near_sgt neb_fem neb_male neuro /// 
	   open phys_act pollen read relig ///
	   risk sat_fam sat_fin sat_frnd sat_job ///
	   self_math swb	   
	   
***********************************************************************************
*** RUN DOFILES
***********************************************************************************

*** Clean Phenotype Data
global RSCRIPT_PATH ""
rscript using "${syntax}/A_phenotypes.R"

*** Clean Phenotype Data
do "${syntax}/B_phenotypes.do"

***Rho Estimation
do "$syntax/C_rho.do"

*** Begin Cleaning WLS Data & Auxillary Files
do "${syntax}/D_clean.do"

*** Descriptive Statistics Table
do "${syntax}/G_tab1.do"

*** FE Beta x PD Beta Scatter (Figure 1)
do "${syntax}/H_fig1.do"

*** PGI on Lifespan (Figure 2)
do "${syntax}/I_fig2.do"

*** SD Test Supplementary Table
do "$syntax/J_tabS2.do"

*** Rho Data to .xlsx
do "$syntax/K_tabS3.do"

*** PGI on Survival to 75 Supplementary Figure
do "$syntax/L_figS1.do"

*** Construct Precision and Bias Figures
global RSCRIPT_PATH ""
rscript using "${syntax}/M_figS23.R"
