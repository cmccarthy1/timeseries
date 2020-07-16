// The following are updated versions of the timeseries utils provided in the original code base,
// some of these changes may be slightly less efficient than the original versions
// but in the cases where there has been a performance hit this has been done for the 
// sake of code scalability/readability/flexibility with AutoML/Analyst and support in mind.

// The following are a number of variable definitions that occur throughout this
// file and are provided at this point to limit repetition 
/* data = dataset on which functionality is to be applied (vector/table/dict etc)
/* lag  = long indicating the number of steps to 'lag' a dataset

\d .tm

.tm.loadfile`:code/optim.q
// loading utilities here for testing puposes (should be wrapped in final version)
\l ml/ml.q
.ml.loadfile`:init.q

// AR/ARMA/SARMA model utilities
/* endog = endogenous variable (time-series)
/* exog  = exogenous variables (additional variables)
/* d     = dictionary containing p,q,tr and seasonal P Q values
/* typ   = typ of model, ARMA or SARMA

// Fit an (S)ARMA model using the Hannan-Rissanen method
/. r       > dictionary of params and data for future predictions
i.SARMAmdl:{[endog;exog;d;typ]
  // Convert exog to matrix
  if[98h=type exog;exog:"f"$i.mat[exog]];
  // Number of lags which must be accounted for
  n:1+max d`p`q;
  // Construct an AR model to estimate the residual error parameters
  errs:i.esterrs[endog;exog;n];
  // Based on the residual errors calculated the coefficients of the ARMA model
  coeff:i.estparam[endog;exog;errs`err;d];
  // if SARIMA model, then use coeffs as starting parameters to calculate sarima coeffs
  if[not all null raze d[`P`Q];coeff:i.SARMA_coeff[endog;exog;errs`err;coeff;d]];
  // Construct the dictionary for return from the function;
  key_vals:get[".tm.i.",typ,"keys"][];
  params:(coeff(::;d[`tr]-1;d[`tr]+til count exog 0)),get[".tm.i.",typ,"params"][endog;coeff;d;errs;n];
  key_vals!params
  }

// Predict future values using an (S)ARMA model
/. r        > future predictions
i.SARMApred:{[mdl;exog;len;typ]
  // allow null to be provided as exogenous variable
  if[exog~(::);exog:()];
  if[not count[mdl`exog_param]~count exog[0];i.err.exog[]];  
  // convert exogenous variable to a matrix if required
  if[98h~type exog;exog:"f"$i.mat exog];
  // if any residual values are present then use these in prediction otherwise use ARpred
  $[count raze mdl[`pred_dict];
    last{x>count y 2}[len;]i.sngpred[mdl`params;exog;mdl`pred_dict;;mdl`estresid;typ]/(mdl`lags;mdl`resid;());
     ARpred[mdl;exog;len]]
  }

// Estimate ARMA model parameters using OLS
/. r       > estimated parameters
i.estparam:{[endog;exog;errors;d]
  // Create lagged matrices for the endogenous variable and residual errors
  endogm:i.lagmat[endog;d`p];
  resid :i.lagmat[errors;d`q];
  // Collect the data needed for estimation
  vals:(exog;endogm;resid);
  // How many data points are required
  m:neg min raze(count[endog]-d[`p`P]),count[errors]-d[`q`Q];
  x:(,'/)m#'vals;
  // add seasonality components
  if[not 0N~d[`P];x:x,'(m #flip[d[`P]xprev\:endog])];
  if[not 0N~d[`Q];x:x,'(m #flip[d[`Q]xprev\:errors])];
  // If required add a trend line variable
  if[d`tr;x:1f,'x];
  y:m#endog;
  first enlist[y]lsq flip x
  }

// Durbin Levinson function to calculate the coefficients in a pure AR model with no trend for a univariate dataset
// Implementation can be found here https://www.stat.purdue.edu/~zhanghao/STAT520/handout/DurbinLevHandout.pdf
/. r - returns coefficients for lag values
i.durbin_lev:{[data;p]
 mat:(1+p;1+p)#0f;
 v:(1+p)#0f;
 mat[1;1]:.tm.i.acf[data;1];
 v[1]:var[data]*(1-xexp[mat[1;1];2]);
 reverse 1_last first(p-1){[data;d]
 mat:d[0];v:d[1];n:d[2];
 k:n+1;
 mat[k;k]:(.tm.i.lagcov[data;k]-sum mat[n;1+til n]mmu .tm.i.lagcov[data]each k-1+til n)%v[n];
 upd:{[data;n;mat;j]
  mat[n;j]-(mat[n+1;n+1]*mat[n;1+n-j])
 }[data;n;mat]each 1+til n;
 mat[k;1+til n]:upd;
 v[k]:v[n]*(1-xexp[mat[k;k];2]);
 (mat;v;n+1)}[data]/(mat;v;1)
 }

// Estimate errors for Hannan Riessanan method
/* n     = AR param to use
/. r     > the residual errors and paramaters used to calculate them
i.esterrs:{[endog;exog;n]
 // Construct an AR model to estimate the residual error parameters
 estresid:ARfit[endog;exog;n;0b]`params;
 // Convert the endogenous variable to lagged matrix
 endogm:i.lagmat[endog;n];
 // Predict future values based on estimations from AR model and use to estimate error
 err:(n _endog)-((neg[count endogm]#exog),'endogm)mmu estresid;
 `params`err!(estresid;err)}

// Use the estimated coefficients as starting points to calculate the sarima coeffs
/. r  > returns the updated coefficients
i.SARMA_coeff:{[endog;exog;resid;coeff;d]
 // data length to use
 len_q:count[resid]-max raze d[`q`Q`seas_add_Q];
 len_p:count[endog]-max raze d[`p`P`seas_add_P];
 // prediction values
 d[`real]:#[m:neg min len_p,len_q;endog];
 // get lagged values
 lag_val:i.lagmat[endog;d`p];
 // get seasonal lag values
 seas_lag:flip d[`P]xprev\:endog;
 // get additional seasonal lag values
 d[`seas_lag_add]:$[(d`p)&min count d`P;#[m;flip d[`seas_add_P]xprev\:endog];2#0f];
 // get resid vals
 resid_val:i.lagmat[resid;d`q];
 seas_resid:flip d[`Q]xprev\:resid;
 d[`seas_resid_add]:$[(d`q)&min count d`Q;#[m;flip d[`seas_add_Q]xprev\:resid];2#0f];
 // normal arima vals
 vals:(exog;lag_val;resid_val;seas_lag;seas_resid);
 d[`norm_mat]:(,'/)m#'vals;
 opt_d:`xk`args!(coeff;d);
 // use optimizer function to improve SARMA coefficients
 optim_q:optimize[i.SARMA_max_likeli;opt_d]`xk}

// prep function for sarima model predictions
/* params = estimated params 
/* d      = dict with p,q,P,Q info
/. r      > returns dictionary with new seasonal additional params
i.prep_SARMA:{[params;d]
  // split up the coefficients to their respective p,q,P,Q parts
  lag_p:(d[`tr] _params)[til d`p];
  lag_q:((d[`tr]+d`p)_params)[til d`q];
  lag_seas_p:((d[`tr]+sum d`q`p)_params)[til count[d`P]];
  lag_seas_q:neg[count d`Q]#params;
  // get the additional seasonal multiplied coefficients
  // These coefficients multiply p x P vals and q x Q vals
  lag_mul:$[(d`p)&min count d`P;(*/)flip lag_p cross lag_seas_p;2#0f];
  resid_mul:$[(d`q)&min count d`Q;(*/)flip  lag_q cross lag_seas_q;2#0f];
  // return dictionary of the additional coefficients
  `add_lag_param`add_resid_param!(lag_mul;resid_mul)}

// function to be passed to maximum likelihood function to calculate SARIMA coefficients
/* params = the parameters for the model
/* d      = dictionary of any additional arguments needed
/. r      > returns the sqrt sum of squared errors
i.SARMA_max_likeli:{[params;d]
 // get additional seasonal parameters 
 d,:i.prep_SARMA[params;d];
 // calculate sarima model including the additional seasonal coeffs
 preds:i.SARMA_pred[params;d];
 // calculate error
 sqrt sum n*n:preds-d`real}

// Predict a single AR/ARMA/SARMA value
/* params   = fit parameters
/* d        = dictionary containing indexing information
/* pvals    = list containing p lag values, q residual errors and previously predicted values
/* estres   = estimated residual error
/* typ      = typ of model bring used
/. r        > list of lag values, preds and predicted values
i.sngpred:{[params;exog;d;pvals;estresid;typ]
  // entry in table from which prediction is being made
  exogval:exog count pvals 2;
  // check if sarma or arma model
  pred:get[".tm.i.sng_",typ,"_pred"][params;pvals;exogval;d];
  // if MA>0, estimate resid errors and append
  if[count pvals 1;estvals:exogval,$["SARMA"~typ;#[neg[d`n];];]pvals[0];
   pvals[1]:(1_pvals[1]),pred-mmu[estresid;estvals]];  
   // append new lag values, for next step calculations
  ((1_pvals[0]),pred;pvals[1];pvals[2],pred)}

// make single prediction for SARIMA model
/* params = parameters for predictions
/* d      = dictionary with appropriate split 
/. r      > returns value for single prediction 
i.sng_SARMA_pred:{[params;pvals;exogval;d]
 // include the additional seasonal parameters
 d,:i.prep_SARMA[params;d];
 d[`seas_resid_add]:$[(d`q)&min count d`Q;pvals[1]d[`seas_add_Q];2#0f];
 d[`seas_lag_add]:$[(d`p)&min count d`P;pvals[0]d[`seas_add_P];2#0f];
 d[`norm_mat]:exogval,i.SARMAval[pvals;d];
 i.SARMA_pred[params;d]}

i.sng_ARMA_pred:{[params;pvals;exogval;d]
 m:exogval,i.ARMAval[pvals;d];
 $[d[`tr];
    // add trend value to offset prediction (OoO not important as 2 single vectors)
    params[0]+m mmu 1_params; 
    // no trend value
    params mmu m]}

// sarima prediction function
/* params = parameters for predictions
/* dict   = dictionary with data appropriately split
/. r       > returns the predicted value
i.SARMA_pred:{[params;d]$[d[`tr];params[0]+;]mmu[d[`norm_mat];d[`tr] _params]+
    mmu[d[`seas_resid_add];d`add_resid_param]+mmu[d[`seas_lag_add];d`add_lag_param]}

// Extract fitted ARMA model params to return
/. r         > list of params needed for future predictions
i.ARMAparams:{[endog;coeff;d;errs;n]
 (d[`p]#neg[sum d`q`p]#coeff;neg[d`q]#coeff),
     (neg[n]#endog;neg[d`q]#errs`err;errs`params),enlist d
 }

// Extract fitted SARIMA model params to return
/. r         > list of params needed for future predictions
i.SARMAparams:{[endog;coeff;d;errs;n]
 // number of seasonal components
 ns:count raze d`P`Q;
 // Separate coeffs into normal and seasonal componants
 coefn:neg[ns]_coeff;coefs:neg[ns]#coeff;
 params:(d[`p]#neg[sum d`q`p]#coefn;neg[d`q]#coefn;count[d`P]#coefs;
   neg count[d`Q]#coefs),((neg n|max[raze d`P`seas_add_P])#endog;
  (neg max raze d`p`Q`seas_add_Q)#errs`err;errs`params);
 // Update dictionary values for seasonality funcs
 d[`P`Q`seas_add_P`seas_add_Q]:d[`P`Q`seas_add_P`seas_add_Q]-min d[`m];
 params,enlist d,`tr`n!d[`tr],n
 }

// Extract appropriate lag and resid values for ARMA future predictions
/. r       > list of lag and resid values in appropriate order for predictions
i.ARMAval:{[pvals;d] raze#[neg[d`p];pvals[0]],pvals[1]}

// Extract appropriate lag and resid values for SARIMA future predictions
/. r      > list of lag and resid values in appropriate order for predictions
i.SARMAval:{[pvals;d];
  raze#[neg[d`p];pvals[0]],#[neg[d`q];pvals[1]],
   pvals[0][d[`P]],pvals[1][d[`Q]]}

// Keys of ARMA and SARIMA model dictionaries
i.ARMAkeys:`params`tr_param`exog_param`p_param`q_param`lags`resid`estresid`pred_dict;
i.SARMAkeys:`params`tr_param`exog_param`p_param`q_param`P_param`Q_param`lags`resid`estresid`pred_dict;

// Akaike Information Criterion

/* true   = true values
/* pred   = predicted values
/* params = set of parameters used in creation of the model
/. r      > aic score
i.aicscore:{[true;pred;params]
  // Calculate residual sum of squares, normalised for number of values
  rss:{wsum[x;x]%y}[true-pred;n:count pred];
  // Number of parameter
  k:sum params;
  aic:(2*k)+n*log rss;
  // if k<40 use the altered aic score
  $[k<40;aic+(2*k*k+1)%n-k-1;aic]
  }

// Fit a model, predict the test, return AIC score for single set of input params
/* train  = Training data as a dictionary with endog and exog data
/* test   = Testing data as a dictionary with endog and exog data
/* len    = Number of steps in the future to be predicted
/* params = parameters used in prediction
/. r      > aic score
i.aicfitscore:{[train;test;len;params]
  // Fit an model using the specified parameters
  mdl :ARIMAfit[train`endog;train`exog;;;;]. params`p`d`q`tr;
  // Predict using the fitted model
  pred:ARIMApred[mdl;test`exog;len];
  // Score the predictions
  i.aicscore[len#test`endog;pred;params]
  }


// Autocorrelation functionality

// Lagged covariance functionality
/. r     > covariance between a dataset at time t and t-lag 
i.lagcov:{[data;lag]cov[neg[lag] _ data;lag _ data]}

// Calculate the autocorrelation between a series and lagged version of itself
i.acf:{[data;lag]i.lagcov[data;lag]%var data}

// Matrix creation/manipulation functionality

// Create a lagged matrix with each row containing 'lag' values
i.lagmat:{[data;lag]data til[count[data]-lag]+\:til lag}

// Matrix from table
i.mat:{[data]flip value flip data}


// Stationarity functionality used to test if datasets are suitable for application of the ARIMA
// and to facilitate transformation of the data to a more suitable form if relevant

// Function to calculate relevant augmented dickey fuller statistics
/* dtype = type of the dataset that's being passed to the function
/. r     > all relevant scores from an augmented dickey fuller test
i.statscores:{[data;dtype]
  // Calculate the augmented dickey-fuller scores for a dict/tab/vector input
  scores:{.ml.fresh.i.adfuller[x]`}@'
    $[98h=dtype;flip data;
      99h=dtype;data;
      dtype in(6h;7h;8h;9h);enlist data;
      '"Inappropriate type provided"];
  flip{x[0 1],(0.05>x 1),value x 4}each$[dtype in 98 99h;value::;]scores
  }

// Are all of the series provided by a user stationary, determined using augmented dickey fuller?
// dickey fuller test
/. r    > boolean indicating if all time series are stationary or not
i.stat:{[data](all/)i.statscores[data;type data][2]}

// time-series differencing to establish stationarity
/* d    = order of time series differencing
/. r    > differenced time series with the first d elements removed
i.diff:{[data;d]d _d{deltas x}/data}

// Seasonal differencing
/*m = order of the seasonal component
/*d = data to apply differencing on
/.r > seasonal diffenced data
i.sdiff:{[m;d][m]_ d-(m xprev d)}

// Revert season differenced data
/*origd  = original data before being differenced
/*dfdata = differenced data
/.r      > the data reverted back to its original format before differencing 
i.revseasdf:{[origs;dfdata]
 seasd:origs,dfdata;
 n:count origs;
 [n]_first{x[1]<y}[;count[seasd]]{[n;sdi]
 sd:sdi[0];i:sdi[1];
 sd[i]:sd[i-n]+sd[i];
 (sd;i+1)}[n]/(seasd;n)}

// Error calls
i.err.steps:{'`$"Exog length not long enough"}
i.err.stat:{'`$"Time series not stationary, try another value of d"}
i.err.len:{'`$"Endog length less than length"}
i.err.exog:{'`$"Test exog length does not match train exog length"}
