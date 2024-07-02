***********************************************************************************
***
***********************************************************************************

***
use "${analytic}", clear

***********************************************************************************
*** PANEL A
***********************************************************************************
preserve

	keep if dna_sib==1

	***make and export summary statistic table
	estimates clear

	***put grads before sibs in table
	replace is_grad = abs(is_grad - 1)

	bysort is_grad: eststo: estpost sum female born dead alive_by_75 imp_span
									   
	esttab /// 
		using "${table}\tab1a_${date}.tex", ///
		tex label nodepvar noobs nonumber plain replace ///
		cells("mean(fmt(2)) sd(fmt(2)) count(fmt(0))") /// 
		collabels("Mean" "SD" "N") ///
		mtitle("Graduate" "Not Graduate")
		
restore	

***********************************************************************************
*** PANEL B
***********************************************************************************

preserve

keep if dna_sib==0

replace dna=1-abs(dna)

***make and export summary statistic table
estimates clear
bysort dna: eststo: estpost sum	is_grad female born dead alive_by_75 imp_span
								
esttab /// 
	using "${table}\tab1b_${date}.tex", ///
	label nodepvar noobs nonumber plain replace ///
	cells("mean(fmt(2)) sd(fmt(2)) count(fmt(0))") /// 
	collabels("Mean" "SD" "N") ///
	mtitle("Genotyped" "Not Genotyped")
restore
