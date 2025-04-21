//Part1
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"
tsmktim time, start(1981m1)
gen t=_n

//ac and pac plots
//ac ld1_cpi_all_can
//pac ld1_cpi_all_can

//1-step ahead
//AR model for 1 step ahead forecast
* Set the desired values for starting and ending t and AR range
local start_t = 33
local end_t = 268
local ar_start = 1
local ar_end = 20

* Create an empty local macro to later store the AR model names for the estimates stats command
local models ""

* Loop through the AR orders and run the regressions
forvalues i = `ar_start'/`ar_end' {
	local j = `i' + 11
    reg ld12_cpi_all_can L(12/`j').ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
    estimates store AR`i'
    * Append the model name to the models macro
    local models "`models' AR`i'"
}

estimates stats `models'


//Best AR based on AIC: AR15, Best AR based on BIC: AR12
//new starting t and ending t to accomodate 18 lags + horizon shift
local start_t = 31
local end_t = 268

//Building ADL models: Granger Causality Test, Based on AIC for AR
disp 0.75*236^(1/3) /*SW lag default for HAC*/

forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').m_base1 if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').m_base1 
} 
//significant at up to 15, we can include up to 15 to test

forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').m2p if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').m2p 
} //significant up till 18, we test up to 18 lags of m2p


forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').m3 if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').m3
} //significant to 18, include 18 lags of m3


forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').bank_rate_l if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').bank_rate_l
} //bank rate l significant up to lag 18, we test till 18


forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').tbill_6mbank_rate if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').tbill_6mbank_rate
} //significant in a few lags including lag 18 we test up to 18 lags of tbill_6mbank_rate


forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').ippi_can  if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').ippi_can 
} //not significant for any lags of ippi_can


forvalues i = 1/18 {
local j = `i'+11 
qui newey ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/`j').unemp_can  if t >= `start_t' & t <= `end_t', lag(4)
testparm L(12/`j').unemp_can 
} //significnat at 18 lags of unemp_can test up to 18 lags


reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR15

reg ld12_cpi_all_can L(12/23).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR12

reg ld12_cpi_all_can L(12/12).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR1

local models_mbase ""

//testing to find best m_base1
local ar_range = "15"
local adl_range = "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+11
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld12_cpi_all_can L(12/`i').ld1_cpi_all_can L(12/`adl').m_base1 if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-11
        estimates store ADL_mbase_`j'
        * Append the model name to the models macro
        local models_mbase "`models_mbase' ADL_mbase_`j'"
    }
}

estimates stats AR15 `models_mbase' // no need to include mbase1 does not improve AR15


local models_m2p ""


//testing to find best m2p
local ar_range = "15"
local adl_range = "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+11
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld12_cpi_all_can L(12/`i').ld1_cpi_all_can L(12/`adl').m2p if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-11
        estimates store ADL_m2p_`j'
        * Append the model name to the models macro
        local models_m2p "`models_m2p' ADL_m2p_`j'"
    }
}

estimates stats AR15 `models_m2p' // include 7 lags of m2p, 12 to 18


local models_br ""


//testing to find best br
local ar_range = "15"
local adl_range = "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+11
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld12_cpi_all_can L(12/`i').ld1_cpi_all_can L(12/`adl').bank_rate_l if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-11
        estimates store ADL_br_`j'
        * Append the model name to the models macro
        local models_br "`models_br' ADL_br_`j'"
    }
}

estimates stats AR15 `models_br' // include 11 lags of br

local models_tbill ""


//testing to find best tbill
local ar_range = "15"
local adl_range = "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+11
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld12_cpi_all_can L(12/`i').ld1_cpi_all_can L(12/`adl').tbill_6mbank_rate if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-11
        estimates store ADL_tbill_`j'
        * Append the model name to the models macro
        local models_tbill "`models_tbill' ADL_tbill_`j'"
    }
}

estimates stats AR15 `models_tbill' // include 18 lags of tbill



local models_unemp ""

//testing to find best unemp
local ar_range = "15"
local adl_range = "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+11
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld12_cpi_all_can L(12/`i').ld1_cpi_all_can L(12/`adl').unemp_can if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-11
        estimates store ADL_unemp_`j'
        * Append the model name to the models macro
        local models_unemp "`models_unemp' ADL_unemp_`j'"
    }
}

estimates stats AR15 `models_unemp' // include 17 lags of unemp




local models_m3 ""

//testing to find best unemp
local ar_range = "15"
local adl_range = "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+11
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld12_cpi_all_can L(12/`i').ld1_cpi_all_can L(12/`adl').m3 if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-11
        estimates store ADL_m3_`j'
        * Append the model name to the models macro
        local models_m3 "`models_m3' ADL_m3_`j'"
    }
}

estimates stats AR15 `models_m3' // include 18 lag of m3


* Define the list of variables with their respective lags
local m2_chosen "L(12/18).m2p"
local br_chosen "L(12/22).bank_rate_l"
local tbill_chosen "L(12/29).tbill_6mbank_rate"
local un_chosen "L(12/28).unemp_can"
local m3_chosen "L(12/29).m3" //L(12/26) default

* Test all combinations of 2 variables
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `m2_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2br
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `tbill_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2tbill
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2un
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2m3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `tbill_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brtbill
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brun
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `tbill_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store tbillun
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `tbill_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store tbillm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store unm3

* Test all combinations of 3 variables
* Test all combinations of 3 variables
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `tbill_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brtbill
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brun
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `tbill_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2tbillun
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `tbill_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2tbillm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2unm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `tbill_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brtbillun
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `tbill_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brtbillm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brunm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `tbill_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store tbillunm3

* Test all combinations of 4 variables
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `tbill_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brtbillun
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `tbill_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brtbillm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brunm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `tbill_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2tbillunm3
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `br_chosen' `tbill_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brtbillunm3

* Test the model with all 5 variables
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can `m2_chosen' `br_chosen' `tbill_chosen' `un_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brtbillunm3

estimates stat AR1 AR12 AR15 ADL_m2p_7 ADL_br_11 ADL_tbill_18 ADL_unemp_17 ADL_m3_18 ///
m2br m2tbill m2un m2m3 brtbill brun brm3 tbillun tbillm3 unm3 ///
m2brtbill m2brun m2brm3 m2tbillun m2tbillm3 m2unm3 brtbillun brtbillm3 brunm3 tbillunm3 ///
m2brtbillun m2brtbillm3 m2brunm3 m2tbillunm3 brtbillunm3 ///
m2brtbillunm3  //also check if same number of observations of all models

//best 3 ADL model is 1) m2br 2) m2un 3) m2brun
estimates stat AR15 m2br m2un m2brun
//Part2
////////////////////////////////////////////////////////////////////////////////
//now we test using the validation set, get ready for forecast combination

//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can ld12_cpi_all_can unemp_can m3 m2p bank_rate_l tbill_6mbank_rate

tsmktim time, start(1981m1)
gen t=_n
local start_t = 31 //restrict start to make all models have same number of observations, for up to AR(15), and accomodating 18 lags.
local end_t = 268
local validation_start = 269
local validation_end = 368
local windowsize = `validation_end' - `validation_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h12 = (L12.ld12_cpi_all_can + L13.ld12_cpi_all_can + L14.ld12_cpi_all_can + L15.ld12_cpi_all_can + L16.ld12_cpi_all_can + L17.ld12_cpi_all_can + L18.ld12_cpi_all_can + L19.ld12_cpi_all_can + L20.ld12_cpi_all_can + L21.ld12_cpi_all_can + L22.ld12_cpi_all_can + L23.ld12_cpi_all_can)/12 if t>= `validation_start' & t <=`validation_end' 


* Define the AR orders
local ar_orders = "12 26" //for both AR(1) and AR(15)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	local ar_order1 = `ar_order' -11
	gen ar`ar_order1'_h12 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i' -12
        
        * Estimate the AR model for the current training window
        qui reg ld12_cpi_all_can L(12/`ar_order').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order1'_h12 = temp_pred if t == `current_end_t' + 1 + 12
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}


//we also generate the forecasts for the models which are the best 3 adl models: 
gen m2br_h12= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -12
	
	* Estimate the AR model for the current training window
	qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/22).bank_rate_l ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2br_h12 = temp_pred if t == `current_end_t' + 1 +12
	
	* Drop the temporary prediction variable
	drop temp_pred
}
//we also generate the forecasts for the models which are the best 3 adl models: 
gen m2un_h12 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -12
	
	* Estimate the AR model for the current training window
	qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/28).unemp_can  ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2un_h12 = temp_pred if t == `current_end_t' + 1 +12
	
	* Drop the temporary prediction variable
	drop temp_pred
}

gen m2brun_h12 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -12
	
	* Estimate the AR model for the current training window
	qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/22).bank_rate_l  L(12/28).unemp_can  ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2brun_h12 = temp_pred if t == `current_end_t' + 1 +12
	
	* Drop the temporary prediction variable
	drop temp_pred
}





//now we do forecast combination - Getting the weights
//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.

//here we just do it for bates-granger and granger ramanathan as they need the true values
//Bates-granger
//first calculate the RMSE
gen true_values = ld12_cpi_all_can if t >= `validation_start' & t <= `validation_end'

egen e_RW_h12 = mean((RW_h12 - true_values)^2)
gen rmse_RW_h12 = sqrt(e_RW_h12)

egen e_ar1_h12 = mean((ar1_h12 - true_values)^2)
gen rmse_ar1_h12 = sqrt(e_ar1_h12)

egen e_ar15_h12 = mean((ar15_h12 - true_values)^2)
gen rmse_ar15_h12 = sqrt(e_ar15_h12)

egen e_m2br_h12 = mean((m2br_h12 - true_values)^2)
gen rmse_m2br_h12 = sqrt(e_m2br_h12)

egen e_m2un_h12 = mean((m2un_h12 - true_values)^2)
gen rmse_m2un_h12= sqrt(e_m2un_h12)

egen e_m2brun_h12 = mean((m2brun_h12 - true_values)^2)
gen rmse_m2brun_h12 = sqrt(e_m2brun_h12)

//Bates-Granger
//calculate variances of the different models to combine
//calculate weights of bates granger in the following way in excel
//1) for every forecast: calculate 1/MSE (not rmse)
//2) then calculate the weight by normalizing by sum across all forecasts
//code to output forecast mse and rmse into excel
scalar mse_RW = e_RW_h12[1]
scalar rmse_RW = rmse_RW_h12[1]

scalar mse_ar1 = e_ar1_h12[1]
scalar rmse_ar1 = rmse_ar1_h12[1]

scalar mse_ar15 = e_ar15_h12[1]
scalar rmse_ar15 = rmse_ar15_h12[1]

scalar mse_m2br= e_m2br_h12[1]
scalar rmse_m2br = rmse_m2br_h12[1]

scalar mse_m2un = e_m2un_h12[1]
scalar rmse_m2un = rmse_m2un_h12[1]

scalar mse_m2brun = e_m2brun_h12[1]
scalar rmse_m2brun= rmse_m2brun_h12[1]

* Create a numeric matrix with MSE and RMSE values
matrix F = (mse_RW, rmse_RW \ ///
		   mse_ar1, rmse_ar1 \ ///
           mse_ar15, rmse_ar15 \ ///
           mse_m2br, rmse_m2br \ /// 
           mse_m2un, rmse_m2un \ /// 
           mse_m2brun, rmse_m2brun) 
* Label the columns
matrix colnames F = mse rmse
* Label the rows with the forecast names
matrix rownames F = RW_h12 ar1_h12 ar15_h12 m2br m2un m2brun

* List the matrix to view it
matrix list F
* Export matrix F to Excel
putexcel set "Bates_Granger_weights_h12.xlsx", modify
putexcel A1=matrix(F), names

//weights for Bates-Granger
//RW = 0.190107828, AR1 = 0.238121124
//AR15 = 0.247281873, m2br = 0.115784164
//m2un = 0.101758675 m2brun = 0.106946336

//Granger-Ramanathan combination
*Unconstrained regression first to see which forecasts are collinear. If find any, drop them:
reg true_values RW_h12 ar1_h12 ar15_h12 m2br_h12 m2un_h12 m2brun_h12, noconstant
*Actually we find none, so we proceed to set the constraint, and use cnsreg to run the regression:
constraint 1 RW_h12 + ar1_h12 + ar15_h12 + m2br_h12 + m2un_h12 + m2brun_h12 = 1
cnsreg true_values RW_h1 ar1_h12 ar15_h12 m2br_h12 m2un_h12 m2brun_h12, constraints(1) noconstant
//dropping negative weights
constraint 2 RW_h12+ar1_h12+ar15_h12 +m2br_h12 +m2brun_h12 = 1
cnsreg true_values RW_h12 ar1_h12 ar15_h12 m2br_h12 m2brun_h12, constraints(2) noconstant
//dropping negative weights
constraint 3 RW_h12 + ar1_h12 + ar15_h12 +m2br_h12 = 1
cnsreg true_values RW_h12 ar1_h12 ar15_h12 m2br_h12, constraints(3) noconstant
constraint 4 RW_h12 + ar1_h12 + ar15_h12 = 1
cnsreg true_values RW_h12 ar1_h12 ar15_h12, constraints(4) noconstant
constraint 5 ar1_h12 + ar15_h12 = 1
cnsreg true_values ar1_h12 ar15_h12, constraints(5) noconstant

//Granger-Ramanathan weights
// ar1_h12 = .3907031, ar15_h12 = .6092969

//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.




//Part 3: same variable names, but now with all data up till before test period
////Final run of everything:
//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can ld12_cpi_all_can unemp_can m2p bank_rate_l 

tsmktim time, start(1981m1)
gen t=_n
local start_t = 31 //restrict start to make all models have same number of observations, for up to AR(15)
local end_t = 368 //have all the data in for the test set
local test_start = 369
local test_end = 468
local windowsize = `test_end' - `test_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast'
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h12 = (L12.ld12_cpi_all_can + L13.ld12_cpi_all_can + L14.ld12_cpi_all_can + L15.ld12_cpi_all_can + L16.ld12_cpi_all_can + L17.ld12_cpi_all_can + L18.ld12_cpi_all_can + L19.ld12_cpi_all_can + L20.ld12_cpi_all_can + L21.ld12_cpi_all_can + L22.ld12_cpi_all_can + L23.ld12_cpi_all_can)/12 if t>= `test_start' & t <=`test_end'  //note: for higher horize


* Define the AR orders
local ar_orders = "12 26" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	local ar_order1 = `ar_order'-11
	gen ar`ar_order1'_h12 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i' - 12
        
        * Estimate the AR model for the current training window
        qui reg ld12_cpi_all_can L(12/`ar_order').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order1'_h12 = temp_pred if t == `current_end_t' + 1 + 12
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}

//we also generate the forecasts for the models which are the best 3 adl models: 
gen m2br_h12= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 12
	
	* Estimate the AR model for the current training window
	qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/22).bank_rate_l ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2br_h12 = temp_pred if t == `current_end_t' + 1 +12
	
	* Drop the temporary prediction variable
	drop temp_pred
}

gen m2un_h12= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 12
	
	* Estimate the AR model for the current training window
	qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/28).unemp_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2un_h12 = temp_pred if t == `current_end_t' + 1 +12
	
	* Drop the temporary prediction variable
	drop temp_pred
}



//we also generate the forecasts for the models which are the best 3 adl models:
gen m2brun_h12 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -12
	
	* Estimate the AR model for the current training window
	qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/22).bank_rate_l L(12/28).unemp_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2brun_h12 = temp_pred if t == `current_end_t' + 1 + 12
	
	* Drop the temporary prediction variable
	drop temp_pred
}

//generate forecast combinations:
//Forecast averaging:
gen forecast_average_h12 = (RW_h12 + ar1_h12 + ar15_h12 + m2br_h12 + m2un_h12 + m2brun_h12)/6

//Bates-Granger Combination:
//using weights from part 2
//weights for Bates-Granger
//RW = 0.190107828, AR1 = 0.238121124
//AR15 = 0.247281873, m2br = 0.115784164
//m2un = 0.101758675 m2brun = 0.106946336
gen bate_granger_h12 = 0.190107828*RW_h12 + 0.238121124*ar1_h12 + 0.247281873*ar15_h12 + 0.115784164*m2br_h12 + 0.101758675*m2un_h12 + 0.106946336*m2brun_h12

//Granger Ramanathan Combination:
//using weights from part 2
// ar1_h12 = .3907031, ar15_h12 = .6092969
gen granger_raman_h12 = 0.3907031*ar1_h12 + 0.6092969*ar15_h12


local start_t = 31 
local end_t = 368 -12 

//BMA of AR and ADL models
* Get the ICs of the different models
qui reg ld12_cpi_all_can L(12/12).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR1
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR15
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/22).bank_rate_l if t >= `start_t' & t <= `end_t', robust
estimates store m2br_h12
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/28).unemp_can if t >= `start_t' & t <= `end_t', robust
estimates store m2un_h12
qui reg ld12_cpi_all_can L(12/26).ld1_cpi_all_can L(12/18).m2p L(12/22).bank_rate_l L(12/28).unemp_can if t >= `start_t' & t <= `end_t', robust
estimates store m2brun_h12

estimates stats AR1 AR15 m2br_h12 m2un_h12 m2brun_h12

//weights following BMA. Copy paste BIC table in stata output into excel and calculate manually (select: right click copy table)
//weights according to BMA:
//AR15: 0.0001362249, m2br_h12: 0.9997963642, m2un_h12:0.0000674109, others are too small, we leave output
gen bma_combined_h12 = 0.0001249455*ar15_h12 + 0.9997963642*m2br_h12 + 0.0000674109*m2un_h12


//WAIC of AR and ADL models
//weights following WAIC. Copy paste AIC table in stata output into excel and calculate manually (select: right click copy table)
// m2br_h12: 0.1449010172, m2un_h12:0.8388659072 m2brun_h12:0.0162330757
// others are too small
gen waic_combined_h12 = 0.1449010172*m2br_h12 + 0.8388659072*m2un_h12 + 0.0162330757*m2brun_h12



//now we do the DM test for each of the 10 forecasts against Random Walk
gen true_values = ld12_cpi_all_can if t >= `test_start' & t <= `test_end'


//errors of forecast
gen e_RW_h12 = true_values - RW_h12 
gen e_ar1_h12 = true_values - ar1_h12
gen e_ar15_h12= true_values - ar15_h12 
gen e_m2br_h12 = true_values - m2br_h12 
gen e_m2un_h12 = true_values - m2un_h12
gen e_m2brun_h12 = true_values - m2brun_h12
gen e_forecast_average_h12 = true_values - forecast_average_h12
gen e_bate_granger_h12 = true_values - bate_granger_h12 
gen e_granger_raman_h12 = true_values - granger_raman_h12 
gen e_bma_combined_h12= true_values - bma_combined_h12
gen e_waic_combined_h12= true_values - waic_combined_h12

//mse and rmse of errors
egen mse_RW_h12 = mean(e_RW_h12^2)
gen rmse_RW_h12 = sqrt(mse_RW_h12)

egen mse_ar1_h12 = mean(e_ar1_h12^2)
gen rmse_ar1_h12 = sqrt(mse_ar1_h12)

egen mse_ar15_h12 = mean(e_ar15_h12^2)
gen rmse_ar15_h12 = sqrt(mse_ar15_h12)

egen mse_m2br_h12 = mean(e_m2br_h12^2)
gen rmse_m2br_h12 = sqrt(mse_m2br_h12)

egen mse_m2un_h12 = mean(e_m2un_h12^2)
gen rmse_m2un_h12 = sqrt(mse_m2un_h12)

egen mse_m2brun_h12 = mean(e_m2brun_h12^2)
gen rmse_m2brun_h12 = sqrt(mse_m2brun_h12)

egen mse_forecast_average_h12= mean(e_forecast_average_h12^2)
gen rmse_forecast_average_h12 = sqrt(mse_forecast_average_h12)

egen mse_bate_granger_h12 = mean(e_bate_granger_h12^2)
gen rmse_bate_granger_h12 = sqrt(mse_bate_granger_h12)

egen mse_granger_raman_h12 = mean(e_granger_raman_h12^2)
gen rmse_granger_raman_h12 = sqrt(mse_granger_raman_h12)

egen mse_bma_combined_h12 = mean(e_bma_combined_h12^2)
gen rmse_bma_combined_h12 = sqrt(mse_bma_combined_h12)

egen mse_waic_combined_h12 = mean(e_waic_combined_h12^2)
gen rmse_waic_combined_h12 = sqrt(mse_waic_combined_h12)


//getting loss differentials and the maxlag to use in DM test
gen d_rw_ar1 =e_RW_h12^2 - e_ar1_h12^2
ac d_rw_ar1
corrgram  d_rw_ar1
//RW-AR1 should use lag 8
gen d_rw_ar15 =e_RW_h12^2 - e_ar15_h12^2
ac d_rw_ar15
corrgram  d_rw_ar15
//RW-AR15 should use lag 9
gen d_rw_adl_1 =e_RW_h12^2 - e_m2br_h12^2
ac d_rw_adl_1
corrgram  d_rw_adl_1
//RW-m2br_h12 should use lag 9
gen d_rw_adl_2 =e_RW_h12^2 - e_m2un_h12^2
ac d_rw_adl_2
corrgram  d_rw_adl_2
//RW-m2un_h12 should use lag 9
gen d_rw_adl_3 =e_RW_h12^2 - e_m2brun_h12^2
ac d_rw_adl_3
corrgram  d_rw_adl_3
//RW-m2brun_h12 should use lag 9
gen d_rw_fa =e_RW_h12^2 - e_forecast_average_h12^2
ac d_rw_fa
corrgram  d_rw_fa
//RW-forecast averaging should use lag 9
gen d_rw_bg =e_RW_h12^2 - e_bate_granger_h12^2
ac d_rw_bg
corrgram  d_rw_bg
//RW-bates granger should use lag 9
gen d_rw_gr =e_RW_h12^2 - e_granger_raman_h12^2
ac d_rw_gr
corrgram  d_rw_gr
//RW-granger ramanathan should use lag 9
gen d_rw_bma =e_RW_h12^2 - e_bma_combined_h12^2
ac d_rw_bma
corrgram  d_rw_bma
//RW-bates granger should use lag 9
gen d_rw_waic =e_RW_h12^2 - e_waic_combined_h12^2
ac d_rw_waic
corrgram  d_rw_waic
//RW-granger ramanathan should use lag 9


//dm test against AR1
gen d_ar1_ar15 =e_ar1_h12^2 - e_ar15_h12^2
ac d_ar1_ar15
corrgram  d_ar1_ar15
//AR1-AR15 should use lag 2
gen d_ar1_adl_1 =e_ar1_h12^2 - e_m2br_h12^2
ac d_ar1_adl_1
corrgram  d_ar1_adl_1
//AR1-m2br_h12 should use lag 29
gen d_ar1_adl_2 =e_ar1_h12^2 - e_m2un_h12^2
ac d_ar1_adl_2
corrgram  d_ar1_adl_2
//AR1-m2un_h12 should use lag 29
gen d_ar1_adl_3 =e_ar1_h12^2 - e_m2brun_h12^2
ac d_ar1_adl_3
corrgram  d_ar1_adl_3
//AR1-m2brun_h12 should use lag 29
gen d_ar1_fa =e_ar1_h12^2 - e_forecast_average_h12^2
ac d_ar1_fa
corrgram  d_ar1_fa
//AR1-forecast averaging should use lag 2
gen d_ar1_bg =e_ar1_h12^2 - e_bate_granger_h12^2
ac d_ar1_bg
corrgram  d_ar1_bg
//AR1-bates granger should use lag 2
gen d_ar1_gr =e_ar1_h12^2 - e_granger_raman_h12^2
ac d_ar1_gr
corrgram  d_ar1_gr
//AR1-granger ramanathan should use lag 12
gen d_ar1_bma =e_ar1_h12^2 - e_bma_combined_h12^2
ac d_ar1_bma
corrgram  d_ar1_bma
//AR1-bates granger should use lag 29
gen d_ar1_waic =e_ar1_h12^2 - e_waic_combined_h12^2
ac d_ar1_waic
corrgram  d_ar1_waic
//AR1-granger ramanathan should use lag 29



//compiling results with dm test:
matrix rmse_matrix = J(11 , 3, 0) //10 rows with 3 columns

//compiling rmse (make sure its in order)
matrix rmse_matrix[1, 1] = rmse_RW_h12[1]
matrix rmse_matrix[2, 1] = rmse_ar1_h12[1]
matrix rmse_matrix[3, 1] = rmse_ar15_h12[1]
matrix rmse_matrix[4, 1] = rmse_m2br_h12[1]
matrix rmse_matrix[5, 1] = rmse_m2un_h12[1]
matrix rmse_matrix[6, 1] = rmse_m2brun_h12[1]
matrix rmse_matrix[7, 1] = rmse_forecast_average_h12[1]
matrix rmse_matrix[8, 1] = rmse_bate_granger_h12[1]
matrix rmse_matrix[9, 1] = rmse_granger_raman_h12[1]
matrix rmse_matrix[10, 1] = rmse_bma_combined_h12[1]
matrix rmse_matrix[11, 1] = rmse_waic_combined_h12[1]

//DM test against RW (make sure its in order):
qui dmariano true_values RW_h1 ar1_h12, crit(mse) kernel(bartlett) maxlag(8)
matrix rmse_matrix[2, 2] = r(p)
qui dmariano true_values RW_h1 ar15_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[3, 2] = r(p)
qui dmariano true_values RW_h1 m2br_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[4, 2] = r(p)
qui dmariano true_values RW_h1 m2un_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[5, 2] = r(p)
qui dmariano true_values RW_h1 m2brun_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[6, 2] = r(p)
qui dmariano true_values RW_h1 forecast_average_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[7, 2] = r(p)
qui dmariano true_values RW_h1 bate_granger_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[8, 2] = r(p)
qui dmariano true_values RW_h1 granger_raman_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[9, 2] = r(p)
qui dmariano true_values RW_h1 bma_combined_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[10, 2] = r(p)
qui dmariano true_values RW_h1 waic_combined_h12, crit(mse) kernel(bartlett) maxlag(9)
matrix rmse_matrix[11, 2] = r(p)

//DM test against AR(1) (make sure its in order):
qui dmariano true_values ar1_h12 ar15_h12, crit(mse) kernel(bartlett) maxlag(2)
matrix rmse_matrix[3, 3] = r(p)
qui dmariano true_values ar1_h12 m2br_h12, crit(mse) kernel(bartlett) maxlag(29)
matrix rmse_matrix[4, 3] = r(p)
qui dmariano true_values ar1_h12 m2un_h12, crit(mse) kernel(bartlett) maxlag(29)
matrix rmse_matrix[5, 3] = r(p)
qui dmariano true_values ar1_h12 m2brun_h12, crit(mse) kernel(bartlett) maxlag(29)
matrix rmse_matrix[6, 3] = r(p)
qui dmariano true_values ar1_h12 forecast_average_h12, crit(mse) kernel(bartlett) maxlag(2)
matrix rmse_matrix[7, 3] = r(p)
qui dmariano true_values ar1_h12 bate_granger_h12, crit(mse) kernel(bartlett) maxlag(2)
matrix rmse_matrix[8, 3] = r(p)
qui dmariano true_values ar1_h12 granger_raman_h12, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[9, 3] = r(p)
qui dmariano true_values ar1_h12 bma_combined_h12, crit(mse) kernel(bartlett) maxlag(29)
matrix rmse_matrix[10, 3] = r(p)
qui dmariano true_values ar1_h12 waic_combined_h12, crit(mse) kernel(bartlett) maxlag(29)
matrix rmse_matrix[11, 3] = r(p)

* Name the columns and rows for clarity 
matrix colnames rmse_matrix = "RMSE" "DM RW" "DM AR(1)"
matrix rownames rmse_matrix =  RW_h12 AR1_h12 AR15_h12 ADL1_h12 ADL2_h12 ADL3_h12 FA_h12 BG_h12 GR_h12 BMA_h12 WAIC_h12

matrix list rmse_matrix

putexcel set "12_step_results.xlsx", replace
putexcel A1=matrix(rmse_matrix), names










