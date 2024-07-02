
***********************************************************************************
*** INITIALIZE BKY (2006) SHARPENED Q-VALUE PROGRAM
***********************************************************************************	
	 
capture program drop bky2006

program define bky2006
	 
quietly sum pval
local totalpvals = r(N)

* Sort the p-values in ascending order and generate a variable that codes each p-value's rank
quietly gen int original_sorting_order = _n
quietly sort pval
quietly gen int rank = _n if pval != .

* Set the initial counter to 1
local qval = 1

* Generate the variable that will contain the 
gen bky06_qval = 1 if pval != .

* Set up a loop that begins by checking which hypotheses are rejected at q = 1.000, then checks which hypotheses are rejected at q = 0.999, then checks which hypotheses are rejected at q = 0.998, etc.  The loop ends by checking which hypotheses are rejected at q = 0.001.

while `qval' > 0 {
	* First Stage
	* Generate the adjusted first stage q level we are testing: q' = q/1+q
	local qval_adj = `qval' / ( 1 + `qval' )
	* Generate value q'*r/M
	gen fdr_temp1 = `qval_adj' * rank / `totalpvals'
	* Generate binary variable checking condition p(r) <= q'*r/M
	gen reject_temp1 = ( fdr_temp1 >= pval ) if pval != .
	* Generate variable containing p-value ranks for all p-values that meet above condition
	gen reject_rank1 = reject_temp1 * rank
	* Record the rank of the largest p-value that meets above condition
	egen total_rejected1 = max( reject_rank1 )

	* Second Stage
	* Generate the second stage q level that accounts for hypotheses rejected in first stage: q_2st = q'*(M/m0)
	local qval_2st = `qval_adj' * ( `totalpvals' / ( `totalpvals' - total_rejected1[1] ) )
	* Generate value q_2st*r/M
	gen fdr_temp2 = `qval_2st' * rank / `totalpvals'
	* Generate binary variable checking condition p(r) <= q_2st*r/M
	gen reject_temp2 = ( fdr_temp2 >= pval ) if pval != .
	* Generate variable containing p-value ranks for all p-values that meet above condition
	gen reject_rank2 = reject_temp2 * rank
	* Record the rank of the largest p-value that meets above condition
	egen total_rejected2 = max( reject_rank2 )

	* A p-value has been rejected at level q if its rank is less than or equal to the rank of the max p-value that meets the above condition
	replace bky06_qval = `qval' if rank <= total_rejected2 & rank != .
	* Reduce q by 0.001 and repeat loop
	drop fdr_temp* reject_temp* reject_rank* total_rejected*
	local qval = `qval' - .001
}

quietly sort original_sorting_order

drop original_sorting_order rank
rename bky06_qval qval

end