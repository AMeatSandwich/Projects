//Part1
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"
tsmktim time, start(1981m1)
gen t=_n

//ac and pac plots
//ac ld3_cpi_all_can
pac ld3_cpi_all_can //lag 17 onwards seem to mostly be insignficant, test until AR(20) to be safe

//3-step ahead
//AR model for 3 step ahead forecast
* Set the desired values for starting and ending t and AR range
local start_t = 23
local end_t = 268
local ar_start = 1
local ar_end = 20

* Create an empty local macro to later store the AR model names for the estimates stats command
local models ""

* Loop through the AR orders and run the regressions
forvalues i = `ar_start'/`ar_end' {
	local j = `i' + 2
    qui reg ld3_cpi_all_can L(3/`j').ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
    estimates store AR`i'
    * Append the model name to the models macro
    local models "`models' AR`i'"
}

estimates stats `models'


//Best AR based on AIC: AR12, Best AR based on BIC: AR3
local start_t = 23
local end_t = 268


//Building ADL models: Granger Causality Test, Based on AIC for AR
disp 0.75*246^(1/3) /*SW lag default for HAC*/

forvalues i = 1/12 {
local j = `i'+2
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').m_base1 if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').m_base1 
} 
//no significant lags, we do not include m_base1 


forvalues i = 1/18 {
local j = `i'+2 
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').m2p if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').m2p 
} //significant up till 18, we test up to 18 lags of m2p


forvalues i = 1/12 {
local j = `i'+2 
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').m3 if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').m3
} //significant at 4 and 5, we can test up to 12 lags of m3


forvalues i = 1/12 {
local j = `i'+2 
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').bank_rate_l if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').bank_rate_l
} //bank rate significant at first lag, we test till 12


forvalues i = 1/12 {
local j = `i'+2 
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').tbill_6mbank_rate if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').tbill_6mbank_rate
} //no significant lags, we do not include tbill_6mbank_rate


forvalues i = 1/12 {
local j = `i'+2 
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').ippi_can  if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').ippi_can 
} //not significant for any lags of ippi_can, we do not include ippi_can


forvalues i = 1/12 {
local j = `i'+2 
qui newey ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/`j').unemp_can  if t >= `start_t' & t <= `end_t', lag(4)
testparm L(3/`j').unemp_can 
} //significant from 5 lags to 9 lags of unemp_can, test up to 12 lags


*begin testing to find the best number of distributed lags
local start_t = 23
local end_t = 268

//testing to find best m2p
local models_m2p ""
local ar_range = "12"
local adl_range = "3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+2
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld3_cpi_all_can L(3/`i').ld1_cpi_all_can L(3/`adl').m2p if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-2
        estimates store ADL_m2p_`j'
        * Append the model name to the models macro
        local models_m2p "`models_m2p' ADL_m2p_`j'"
    }
}

estimates stats AR12 `models_m2p' // include 5 lags of m2p, L(3/7)


//testing to find best m3
local models_m3 ""
local ar_range = "12"
local adl_range = "3 4 5 6 7 8 9 10 11 12 13 14"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+2
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld3_cpi_all_can L(3/`i').ld1_cpi_all_can L(3/`adl').m3 if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-2
        estimates store ADL_m3_`j'
        * Append the model name to the models macro
        local models_m3 "`models_m3' ADL_m3_`j'"
    }
}

estimates stats AR12 `models_m3' // include 6 lags of m3, L(3/8)


//testing to find best br
local models_br ""
local ar_range = "12"
local adl_range = "3 4 5 6 7 8 9 10 11 12 13 14"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+2
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld3_cpi_all_can L(3/`i').ld1_cpi_all_can L(3/`adl').bank_rate_l if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-2
        estimates store ADL_br_`j'
        * Append the model name to the models macro
        local models_br "`models_br' ADL_br_`j'"
    }
}

estimates stats AR12 `models_br' // include 1 lag of br, L(3/3)


//testing to find best unemp
local models_unemp ""
local ar_range = "12"
local adl_range = "3 4 5 6 7 8 9 10 11 12 13 14"
* Loop through the orders of AR_ADL
foreach ar in `ar_range' {
	local i = `ar'+2
	
    * Loop through the orders
    foreach adl in `adl_range' {
        qui reg ld3_cpi_all_can L(3/`i').ld1_cpi_all_can L(3/`adl').unemp_can if t >= `start_t' & t <= `end_t', robust
		local j = `adl'-2
        estimates store ADL_unemp_`j'
        * Append the model name to the models macro
        local models_unemp "`models_unemp' ADL_unemp_`j'"
    }
}

estimates stats AR12 `models_unemp' // include 7 lags of unemp, L(3/9)


* Define the list of variables with their respective lags
local m2_chosen "L(3/7).m2p"
local m3_chosen "L(3/8).m3"
local br_chosen "L(3/3).bank_rate_l"
local un_chosen "L(3/9).unemp_can"

local start_t = 23
local end_t = 268

* Test all combinations of 2 variables
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `m3_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2m3
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `br_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2br
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2un
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m3_chosen' `br_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m3br
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m3_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m3un
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store brun

* Test all combinations of 3 variables
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `m3_chosen' `br_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2m3br
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `m3_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2m3un
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2brun
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m3_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m3brun

* Test the combination of 4 variables
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can `m2_chosen' `m3_chosen' `br_chosen' `un_chosen' if t >= `start_t' & t <= `end_t', robust
estimates store m2m3brun


estimates stat AR1 AR3 AR12 ADL_m2p_7 ADL_m2p_5 ADL_m3_6 ADL_br_1 ADL_unemp_7 ///
m2m3 m2br m2un m3br m3un brun ///
m2m3br m2m3un m2brun m3brun ///
m2m3brun //also check if same number of observations of all models

//best 3 ADL model is 1) m2brun 2) m2un 3) m2br
estimates stat AR12 m2brun m2un m2br
//Part2
////////////////////////////////////////////////////////////////////////////////
//now we test using the validation set, get ready for forecast combination

//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can ld3_cpi_all_can m2p bank_rate_l unemp_can 

tsmktim time, start(1981m1)
gen t=_n
local start_t = 23 //restrict start to make all models have same number of observations, for up to AR(20) as tested in part 1 AIC
local end_t = 268
local validation_start = 269
local validation_end = 368
local windowsize = `validation_end' - `validation_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h3 = (L3.ld3_cpi_all_can + L4.ld3_cpi_all_can + L5.ld3_cpi_all_can)/3 if t>= `validation_start' & t <=`validation_end' 


* Define the AR orders
local ar_orders = "3 14" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	local ar_order1 = `ar_order' -2
	gen ar`ar_order1'_h3 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i' -3
        
        * Estimate the AR model for the current training window
        qui reg ld3_cpi_all_can L(3/`ar_order').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order1'_h3 = temp_pred if t == `current_end_t' + 1 + 3
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}


//we also generate the forecasts for the models which are the best 3 adl models: 1st model
gen m2brun_h3= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -3
	
	* Estimate the AR model for the current training window
	qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/3).bank_rate_l L(3/9).unemp_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2brun_h3 = temp_pred if t == `current_end_t' + 1 +3
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 2nd model
gen m2un_h3= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 3
	
	* Estimate the AR model for the current training window
	qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/9).unemp_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2un_h3 = temp_pred if t == `current_end_t' + 1 + 3
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 3rd and last model
gen m2br_h3 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -3
	
	* Estimate the AR model for the current training window
	qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/3).bank_rate_l  ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2br_h3 = temp_pred if t == `current_end_t' + 1 + 3
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//now we do forecast combination - Getting the weights
//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.
local start_t = 23
local end_t = 268
local validation_start = 269
local validation_end = 368
local windowsize = `validation_end' - `validation_start'
//here we just do it for bates-granger and granger ramanathan as they need the true values
//Bates-granger
//first calculate the RMSE
gen true_values = ld3_cpi_all_can if t >= `validation_start' & t <= `validation_end'

egen e_RW_h3 = mean((RW_h3- true_values)^2)
gen rmse_RW_h3 = sqrt(e_RW_h3)

egen e_ar1_h3 = mean((ar1_h3- true_values)^2)
gen rmse_ar1_h3 = sqrt(e_ar1_h3)

egen e_ar12_h3 = mean((ar12_h3- true_values)^2)
gen rmse_ar12_h3 = sqrt(e_ar12_h3)

egen e_m2brun_h3 = mean((m2brun_h3- true_values)^2)
gen rmse_m2brun_h3 = sqrt(e_m2brun_h3)

egen e_m2un_h3 = mean((m2un_h3- true_values)^2)
gen rmse_m2un_h3= sqrt(e_m2un_h3)

egen e_m2br_h3 = mean((m2br_h3- true_values)^2)
gen rmse_m2br_h3 = sqrt(e_m2br_h3)

//Bates-Granger
//calculate variances of the different models to combine
//calculate weights of bates granger in the following way in excel
//1) for every forecast: calculate 1/MSE (not rmse)
//2) then calculate the weight by normalizing by sum across all forecasts
//code to output forecast mse and rmse into excel
scalar mse_RW = e_RW_h3[1]
scalar rmse_RW = rmse_RW_h3[1]

scalar mse_ar1 = e_ar1_h3[1]
scalar rmse_ar1 = rmse_ar1_h3[1]

scalar mse_ar12 = e_ar12_h3[1]
scalar rmse_ar12 = rmse_ar12_h3[1]

scalar mse_m2brun= e_m2brun_h3[1]
scalar rmse_m2brun = rmse_m2brun_h3[1]

scalar mse_m2un = e_m2un_h3[1]
scalar rmse_m2un = rmse_m2un_h3[1]

scalar mse_m2br = e_m2br_h3[1]
scalar rmse_m2br= rmse_m2br_h3[1]

* Create a numeric matrix with MSE and RMSE values
matrix F = (mse_RW, rmse_RW \ ///
		   mse_ar1, rmse_ar1 \ ///
           mse_ar12, rmse_ar12 \ ///
           mse_m2brun, rmse_m2brun \ /// 
           mse_m2un, rmse_m2un \ /// 
           mse_m2br, rmse_m2br) 
* Label the columns
matrix colnames F = mse rmse
* Label the rows with the forecast names
matrix rownames F = RW_h3 ar1_h3 ar12_h3 m2brun m2un m2br

* List the matrix to view it
matrix list F
* Export matrix F to Excel
putexcel set "Bates_Granger_weights_h3.xlsx", modify
putexcel A1=matrix(F), names

//weights for Bates-Granger
//RW = 0.08343851, AR1 = 0.198328899
//AR12 = 0.209297952, m2brun = 0.161326156
//m2un = 0.162304732, m2br = 0.179026919


//Granger-Ramanathan combination
*Unconstrained regression first to see which forecasts are collinear. If find any, drop them:
reg true_values RW_h3 ar1_h3 ar12_h3 m2brun_h3 m2un_h3 m2br_h3, noconstant
*Actually we find none, so we proceed to set the constraint, and use cnsreg to run the regression:
constraint 1 RW_h3 + ar1_h3 + ar12_h3 + m2brun_h3 + m2un_h3 + m2br_h3 = 1
cnsreg true_values RW_h3 ar1_h3 ar12_h3 m2brun_h3 m2un_h3 m2br_h3, constraints(1) noconstant
//dropping negative weights --> drop RW_h3 and m2brun_h3
constraint 2 ar1_h3 + ar12_h3 + m2un_h3 + m2br_h3 = 1
cnsreg true_values ar1_h3 ar12_h3 m2un_h3 m2br_h3, constraints(2) noconstant
//dropping negative weights --> drop m2un_h3
constraint 3 ar1_h3 + ar12_h3 + m2br_h3 = 1
cnsreg true_values ar1_h3 ar12_h3 m2br_h3, constraints(3) noconstant

//Granger-Ramanathan weights
// ar1_h3 = .3300553   , ar12_h3 = 0.6184087  , m2br_h3 = 0.051536

//We do forecast averaging, BMA and WAIC in the final step, because it doesnt need the real values.




//Part 3: same variable names, but now with all data up till before test period
////Final run of everything:
//1-step ahead rolling window for AR models
clear all
cd "`c(pwd)'"
import delimited "edited_balanced_can_md.csv"

keep date ld1_cpi_all_can ld3_cpi_all_can m2p bank_rate_l unemp_can

tsmktim time, start(1981m1)
gen t=_n

local start_t = 23 //restrict start to make all models have same number of observations, for up to AR(20) as tested in part 1 AIC
local end_t = 368 //have all the data in for the test set
local test_start = 369
local test_end = 468
local windowsize = `test_end' - `test_start' //this is one less than the actual window size, to help with looping


//Generate random walk forecast'
//note: for a given, rw is average of past h periods, so for 2 step ahead, forecast of t+3, is average of t, t-1, t-2
//please RW forecast accordingly
gen RW_h3 = (L3.ld3_cpi_all_can + L4.ld3_cpi_all_can + L5.ld3_cpi_all_can)/3 if t>= `test_start' & t <=`test_end'   //note: for higher horize


* Define the AR orders
local ar_orders = "3 14" //for both AR(1) and AR(12)
* Loop over each AR order
foreach ar_order in `ar_orders' {
    * Generate predictions for each AR order using a rolling window
	local ar_order1 = `ar_order'-2
	gen ar`ar_order1'_h3 = .
    forvalues i = 0/`windowsize' {
        * Adjust the start and end for the training window
        local current_start_t = `start_t' + `i'
        local current_end_t = `end_t' + `i' -3
        
        * Estimate the AR model for the current training window
        qui reg ld3_cpi_all_can L(3/`ar_order').ld1_cpi_all_can if t >= `current_start_t' & t <= `current_end_t', robust
        
        * Predict the next value
        predict temp_pred, xb
        
        * Store the forecasted value in the pred_ar variable
        replace ar`ar_order1'_h3 = temp_pred if t == `current_end_t' + 1 + 3
        
        * Drop the temporary prediction variable
        drop temp_pred
    }
}


//we also generate the forecasts for the models which are the best 3 adl models: 1st model
gen m2brun_h3= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' -3 
	
	* Estimate the AR model for the current training window
	qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/3).bank_rate_l L(3/9).unemp_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2brun_h3 = temp_pred if t == `current_end_t' + 1 +3
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 2nd model 
gen m2un_h3= .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 3
	
	* Estimate the AR model for the current training window
	qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/9).unemp_can ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2un_h3 = temp_pred if t == `current_end_t' + 1 + 3
	
	* Drop the temporary prediction variable
	drop temp_pred
}


//we also generate the forecasts for the models which are the best 3 adl models: 3rd and last model
gen m2br_h3 = .
forvalues i = 0/`windowsize' {
	* Adjust the start and end for the training window
	local current_start_t = `start_t' + `i'
	local current_end_t = `end_t' + `i' - 3
	
	* Estimate the AR model for the current training window
	qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/3).bank_rate_l ///
	if t >= `current_start_t' & t <= `current_end_t', robust
	
	* Predict the next value
	predict temp_pred, xb
	
	* Store the forecasted value in the pred_ar variable
	replace m2br_h3 = temp_pred if t == `current_end_t' + 1 + 3
	
	* Drop the temporary prediction variable
	drop temp_pred
}

//generate forecast combinations:
//Forecast averaging:
gen forecast_average_h3 = (RW_h3 + ar1_h3 + ar12_h3 + m2brun_h3 + m2un_h3 + m2br_h3)/6

//Bates-Granger Combination:
//using weights from part 2
//weights for Bates-Granger
//RW = 0.08343851, AR1 = 0.198328899
//AR12 = 0.209297952, m2brun = 0.161326156
//m2un = 0.162304732, m2br = 0.179026919
gen bate_granger_h3 = 0.08343851*RW_h3 + 0.198328899*ar1_h3 + 0.209297952*ar12_h3 + 0.161326156*m2brun_h3 + 0.162304732*m2un_h3 + 0.179026919*m2br_h3

//Granger Ramanathan Combination:
//using weights from part 2
//Granger-Ramanathan weights
// ar1_h3 = .3300553   , ar12_h3 = 0.6184087  , m2br_h3 = 0.051536
gen granger_raman_h3 = 0.3300553*ar1_h3 + 0.6184087*ar12_h3 + 0.051536*m2br_h3


//BMA of AR and ADL models
* Get the ICs of the different models
local start_t = 23 
local end_t = 368 - 3 

qui reg ld3_cpi_all_can L(3/3).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR1
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can if t >= `start_t' & t <= `end_t', robust
estimates store AR12
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/3).bank_rate_l L(3/9).unemp_can if t >= `start_t' & t <= `end_t', robust
estimates store m2brun_h3
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/9).unemp_can if t >= `start_t' & t <= `end_t', robust
estimates store m2un_h3
qui reg ld3_cpi_all_can L(3/14).ld1_cpi_all_can L(3/7).m2p L(3/3).bank_rate_l if t >= `start_t' & t <= `end_t', robust
estimates store m2br_h3

estimates stats AR1 AR12 m2brun_h3 m2un_h3 m2br_h3

//weights following BMA. Copy paste BIC table in stata output into excel and calculate manually (select: right click copy table)
//weights according to BMA:
gen bma_combined_h3 = 0.868081436*ar1_h3 + 0.1314709*ar12_h3 + 0.000447663*m2br_h3


//WAIC of AR and ADL models
//weights following WAIC. Copy paste AIC table in stata output into excel and calculate manually (select: right click copy table)
//others are too small, we leave them out
gen waic_combined_h3 = 0.002250931*ar12_h3 + 0.083000279*m2brun_h3 + 0.148241687*m2un_h3 + 0.766507103*m2br_h3


//now we do the DM test for each of the 10 forecasts against Random Walk
local test_start = 369
local test_end = 468

gen true_values = ld3_cpi_all_can if t >= `test_start' & t <= `test_end'

//errors of forecast
gen e_RW_h3 = true_values - RW_h3
gen e_ar1_h3 = true_values - ar1_h3
gen e_ar12_h3= true_values - ar12_h3 
gen e_m2brun_h3 = true_values - m2brun_h3 
gen e_m2un_h3 = true_values - m2un_h3
gen e_m2br_h3 = true_values - m2br_h3
gen e_forecast_average_h3 = true_values - forecast_average_h3
gen e_bate_granger_h3 = true_values - bate_granger_h3
gen e_granger_raman_h3 = true_values - granger_raman_h3 
gen e_bma_combined_h3= true_values - bma_combined_h3
gen e_waic_combined_h3= true_values - waic_combined_h3

//mse and rmse of errors
egen mse_RW_h3 = mean(e_RW_h3^2)
gen rmse_RW_h3 = sqrt(mse_RW_h3)

egen mse_ar1_h3 = mean(ar1_h3^2)
gen rmse_ar1_h3 = sqrt(mse_ar1_h3)

egen mse_ar12_h3 = mean(ar12_h3^2)
gen rmse_ar12_h3 = sqrt(mse_ar12_h3)

egen mse_m2brun_h3 = mean(m2brun_h3^2)
gen rmse_m2brun_h3 = sqrt(mse_m2brun_h3)

egen mse_m2un_h3 = mean(m2un_h3^2)
gen rmse_m2un_h3 = sqrt(mse_m2un_h3)

egen mse_m2br_h3 = mean(m2br_h3^2)
gen rmse_m2br_h3 = sqrt(mse_m2br_h3)

egen mse_forecast_average_h3= mean(forecast_average_h3^2)
gen rmse_forecast_average_h3 = sqrt(mse_forecast_average_h3)

egen mse_bate_granger_h3 = mean(bate_granger_h3^2)
gen rmse_bate_granger_h3 = sqrt(mse_bate_granger_h3)

egen mse_granger_raman_h3 = mean(granger_raman_h3^2)
gen rmse_granger_raman_h3 = sqrt(mse_granger_raman_h3)

egen mse_bma_combined_h3 = mean(bma_combined_h3^2)
gen rmse_bma_combined_h3 = sqrt(mse_bma_combined_h3)

egen mse_waic_combined_h3 = mean(waic_combined_h3^2)
gen rmse_waic_combined_h3 = sqrt(mse_waic_combined_h3)


//getting loss differentials and the maxlag to use in DM test
gen d_rw_ar1 =e_RW_h3^2 - e_ar1_h3^2
ac d_rw_ar1
corrgram  d_rw_ar1
//RW-AR1 should use lag 13
gen d_rw_ar12 =e_RW_h3^2  - e_ar12_h3^2
ac d_rw_ar12
corrgram  d_rw_ar12
//RW-AR12 should use lag 13
gen d_rw_adl_1 =e_RW_h3^2  - e_m2brun_h3^2
ac d_rw_adl_1
corrgram  d_rw_adl_1
//RW-m2brun_h3 should use lag 13
gen d_rw_adl_2 =e_RW_h3^2  - e_m2un_h3^2
ac d_rw_adl_2
corrgram  d_rw_adl_2
//RW-m2un_h3 should use lag 13,
gen d_rw_adl_3 =e_RW_h3^2  - e_m2br_h3^2
ac d_rw_adl_3
corrgram  d_rw_adl_3
//RW-m2br_h3 should use lag 13
gen d_rw_fa =e_RW_h3^2  - e_forecast_average_h3^2
ac d_rw_fa
corrgram  d_rw_fa
//RW-forecast averaging should use lag 13
gen d_rw_bg =e_RW_h3^2  - e_bate_granger_h3^2
ac d_rw_bg
corrgram  d_rw_bg
//RW-bates granger should use lag 13
gen d_rw_gr =e_RW_h3^2  - e_granger_raman_h3^2
ac d_rw_gr
corrgram  d_rw_gr
//RW-granger ramanathan should use lag 13
gen d_rw_bma =e_RW_h3^2  - e_bma_combined_h3^2
ac d_rw_bma
corrgram  d_rw_bma
//RW- BMA should use lag 13
gen d_rw_waic =e_RW_h3^2  - e_waic_combined_h3^2
ac d_rw_waic
corrgram  d_rw_waic
//RW-WAIC should use lag 13


//dm test against AR1
gen d_ar1_ar12 =e_ar1_h3^2 - e_ar12_h3^2
ac d_ar1_ar12
corrgram  d_ar1_ar12
//AR1-AR12 should use lag 24
gen d_ar1_adl_1 =e_ar1_h3^2 - e_m2brun_h3^2
ac d_ar1_adl_1
corrgram  d_ar1_adl_1
//AR1-m2brun_h3 should use lag 12
gen d_ar1_adl_2 =e_ar1_h3^2 - e_m2un_h3^2
ac d_ar1_adl_2
corrgram  d_ar1_adl_2
//AR1-m2un_h3 should use lag 12
gen d_ar1_adl_3 =e_ar1_h3^2 - e_m2br_h3^2
ac d_ar1_adl_3
corrgram  d_ar1_adl_3
//AR1-m2br_h3 should use lag 12
gen d_ar1_fa =e_ar1_h3^2 - e_forecast_average_h3^2
ac d_ar1_fa
corrgram  d_ar1_fa
//AR1-forecast averaging should use lag 24
gen d_ar1_bg =e_ar1_h3^2 - e_bate_granger_h3^2
ac d_ar1_bg
corrgram  d_ar1_bg
//AR1-bates granger should use lag 24
gen d_ar1_gr =e_ar1_h3^2 - e_granger_raman_h3^2
ac d_ar1_gr
corrgram  d_ar1_gr
//AR1-granger ramanathan should use lag 24
gen d_ar1_bma =e_ar1_h3^2 - e_bma_combined_h3^2
ac d_ar1_bma
corrgram  d_ar1_bma
//AR1-bma should use lag 24
gen d_ar1_waic =e_ar1_h3^2 - e_waic_combined_h3^2
ac d_ar1_waic
corrgram  d_ar1_waic
//AR1-waic should use lag 12



//compiling results with dm test:
matrix rmse_matrix = J(11 , 3, 0) //11 rows with 3 columns

//compiling rmse (make sure its in order)
matrix rmse_matrix[1, 1] = rmse_RW_h3[1]
matrix rmse_matrix[2, 1] = rmse_ar1_h3[1]
matrix rmse_matrix[3, 1] = rmse_ar12_h3[1]
matrix rmse_matrix[4, 1] = rmse_m2brun_h3[1]
matrix rmse_matrix[5, 1] = rmse_m2un_h3[1]
matrix rmse_matrix[6, 1] = rmse_m2br_h3[1]
matrix rmse_matrix[7, 1] = rmse_forecast_average_h3[1]
matrix rmse_matrix[8, 1] = rmse_bate_granger_h3[1]
matrix rmse_matrix[9, 1] = rmse_granger_raman_h3[1]
matrix rmse_matrix[10, 1] = rmse_bma_combined_h3[1]
matrix rmse_matrix[11, 1] = rmse_waic_combined_h3[1]

//DM test against RW (make sure its in order):
qui dmariano true_values RW_h3 ar1_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[2, 2] = r(p)
qui dmariano true_values RW_h3 ar12_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[3, 2] = r(p)
qui dmariano true_values RW_h3 m2brun_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[4, 2] = r(p)
qui dmariano true_values RW_h3 m2un_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[5, 2] = r(p)
qui dmariano true_values RW_h3 m2br_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[6, 2] = r(p)
qui dmariano true_values RW_h3 forecast_average_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[7, 2] = r(p)
qui dmariano true_values RW_h3 bate_granger_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[8, 2] = r(p)
qui dmariano true_values RW_h3 granger_raman_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[9, 2] = r(p)
qui dmariano true_values RW_h3 bma_combined_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[10, 2] = r(p)
qui dmariano true_values RW_h3 waic_combined_h3, crit(mse) kernel(bartlett) maxlag(13)
matrix rmse_matrix[11, 2] = r(p)

//DM test against AR(1) (make sure its in order):
qui dmariano true_values ar1_h3 ar12_h3, crit(mse) kernel(bartlett) maxlag(24)
matrix rmse_matrix[3, 3] = r(p)
qui dmariano true_values ar1_h3 m2br_h3, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[4, 3] = r(p)
qui dmariano true_values ar1_h3 m2un_h3, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[5, 3] = r(p)
qui dmariano true_values ar1_h3 m2br_h3, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[6, 3] = r(p)
qui dmariano true_values ar1_h3 forecast_average_h3, crit(mse) kernel(bartlett) maxlag(24)
matrix rmse_matrix[7, 3] = r(p)
qui dmariano true_values ar1_h3 bate_granger_h3, crit(mse) kernel(bartlett) maxlag(24)
matrix rmse_matrix[8, 3] = r(p)
qui dmariano true_values ar1_h3 granger_raman_h3, crit(mse) kernel(bartlett) maxlag(24)
matrix rmse_matrix[9, 3] = r(p)
qui dmariano true_values ar1_h3 bma_combined_h3, crit(mse) kernel(bartlett) maxlag(24)
matrix rmse_matrix[10, 3] = r(p)
qui dmariano true_values ar1_h3 waic_combined_h3, crit(mse) kernel(bartlett) maxlag(12)
matrix rmse_matrix[11, 3] = r(p)

* Name the columns and rows for clarity 
matrix colnames rmse_matrix = "RMSE" "DM RW" "DM AR(1)"
matrix rownames rmse_matrix =  RW_h3 AR1_h3 AR12_h3 ADL1_h3 ADL2_h3 ADL3_h3 FA_h3 BG_h3 GR_h3 BMA_h3 WAIC_h3

matrix list rmse_matrix

putexcel set "3_step_results.xlsx", replace
putexcel A1=matrix(rmse_matrix), names










