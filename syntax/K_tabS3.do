
do "${syntax}/programs/bky2006.do"

***********************************************************************************
***EXPORT ESTIMATED RHO FOR THE APPENDIX
***********************************************************************************

***
use "${analytic}", clear

keep *rho*

duplicates drop

gen id = 1

reshape long se_rho_pgi_ rho_pgi_, i(id) j(pheno) string

drop id

rename *_ *

***generate pval
gen pval = (2 * ttail(2105), abs(((rho_pgi) - .5)/se_rho_pgi))

format rho_pgi se_rho_pgi pval %04.3f	

***generate qval
drop if pheno=="meta"
bky2006	

order pheno rho se_rho pval qval

***label columns
label var pheno "Polygenic Score"
label var rho_pgi "ρ(g1j,g2j)"
label var se_rho_pgi "Standard Error"
label var pval "p-Value"
label var qval "q-Value (including meta-PGI)"

***label rows
replace pheno = "Adventurousness" if pheno=="adv"
replace pheno = "Age at First Birth" if pheno=="birth"
replace pheno = "Age at First Menses" if pheno=="menses"
replace pheno = "Age Voice Deepened" if pheno=="deep"
replace pheno = "Attention Deficit Hyperactivity Disorder" if pheno=="adhd"
replace pheno = "Allergy: Cat" if pheno=="cat"
replace pheno = "Allergy: Dust" if pheno=="dust"
replace pheno = "Allergy: Pollen" if pheno=="pollen"
replace pheno = "Asthma/Eczema/Rhinitis" if pheno=="aer"
replace pheno = "Asthma" if pheno=="asthma" 
replace pheno = "Alcohol Misuse" if pheno=="alcoh"
replace pheno = "Body Mass Index" if pheno=="bmi"
replace pheno = "Cannabis Use" if pheno=="canna"
replace pheno = "Cognitive Empathy" if pheno=="cog_emp"
replace pheno = "Childhood Reading" if pheno=="read" 
replace pheno = "Chronic Obstructive Pulmonary Disease" if pheno=="copd"
replace pheno = "Cigarettes per Day" if pheno=="cig_day"
replace pheno = "Cognitive Performance" if pheno=="cog"
replace pheno = "Delay Discounting" if pheno=="dly_disc"
replace pheno = "Depressive Symptoms" if pheno=="dep"
replace pheno = "Drinks per Week" if pheno=="drinks"
replace pheno = "Educational Attainment" if pheno=="edu"
replace pheno = "Ever Smoker" if pheno=="ever_smk" 
replace pheno = "Extraversion" if pheno=="extra" 
replace pheno = "Life Satisfaction: Family" if pheno=="sat_fam"
replace pheno = "Life Satisfaction: Finance" if pheno=="sat_fin"
replace pheno = "Life Satisfaction: Friend" if pheno=="sat_frnd" 
replace pheno = "Life Satisfaction: Work" if pheno=="sat_job"
replace pheno = "Hayfever" if pheno=="hay" 
replace pheno = "Height" if pheno=="hgt" 
replace pheno = "Highest Math" if pheno=="high_math" 
replace pheno = "Left Out of Social Activity" if pheno=="leftout"
replace pheno = "Loneliness" if pheno=="lonely"
replace pheno = "Morning Person" if pheno=="chrono"
replace pheno = "Migraine" if pheno=="migrn" 
replace pheno = "Narcissism" if pheno=="narci"
replace pheno = "Nearsightedness" if pheno=="near_sgt" 
replace pheno = "Number Ever Born: Men" if pheno=="neb_male" 
replace pheno = "Number Ever Born: Women" if pheno=="neb_fem" 
replace pheno = "Neuroticism" if pheno=="neuro" 
replace pheno = "Openness" if pheno=="open" 
replace pheno = "Physical Activity" if pheno=="phys_act"
replace pheno = "Religious Attendance" if pheno=="relig"
replace pheno = "Risk Tolerance" if pheno=="risk"
replace pheno = "Self-Rated Health" if pheno=="health"
replace pheno = "Self-Rated Math Ability" if pheno=="self_math"
replace pheno = "Subjective Well-Being" if pheno=="swb"

export excel using "${table}/tabS3_${date}.xlsx", firstrow(varlabels) replace
