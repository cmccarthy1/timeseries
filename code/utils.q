// The following are updated versions of the timeseries utils provided in the original code base,
// some of these changes may be slightly less efficient than the original versions
// but in the cases where there has been a performance hit this has been done for the 
// sake of code scalability/readability/flexibility with AutoML/Analyst and support in mind.

// The following are a number of variable definitions that occur throughout this
// file and are provided at this point to limit repetition 
/* data = dataset on which functionality is to be applied (vector/table/dict etc)
/* lag  = long indicating the number of steps to 'lag' a dataset

\d .tm

// loading utilities here for testing puposes (should be wrapped in final version)
\l ml/ml.q
.ml.loadfile`:init.q

// AR/ARMA model utilities
/* endog = endogenous variable (time-series)
/* exog  = exogenous variables (additional variables)
/* p     = number of autoregressive terms
/* q     = number of lagged forecast errors
/* tr    = is a trend line required

// Fit an (S)ARMA model using the Hannan-Rissanen method
/. r       > dictionary of params and data for future predictions
i.SARMAmdl:{[endog;exog;dict;tr;typ]
  // Convert exog to matrix
  if[98h=type exog;exog:"f"$i.mat[exog]];
  // Number of lags which must be accounted for
  n:1+max dict`p`q;
  // Construct an AR model to estimate the residual error parameters
  errs:i.esterrs[endog;exog;n];
  // Based on the residual errors calculated the coefficients of the ARMA model
  coeff:i.estparam[endog;exog;errs`err;dict;tr];
  // Construct the dictionary for return from the function;
  key_vals:get[".tm.i.",typ,"keys"][];
  params:(coeff(::;tr-1;tr+til count exog 0)),get[".tm.i.",typ,"params"][endog;coeff;dict;errs;n];
  key_vals!params
  }

// Predict future values using an (S)ARMA model
/. r        > future predictions
i.SARMApred:{[mdl;exog;len;typ]
  // allow null to be provided as exogenous variable
  if[exog~(::);exog:()];
  // convert exogenous variable to a matrix if required
  if[98h~type exog;exog:"f"$i.mat exog];
  // if any residual values are present then use these in prediction otherwise use ARpred
  $[count raze mdl[`pred_dict];
    last{x>count y 2}[len;]i.sngpred[mdl`params;exog;mdl`pred_dict;;mdl`estresid;typ]/(mdl`lags;mdl`resid;());
     ARpred[mdl;exog;len]]
  }

// Estimate ARMA model parameters using OLS
/* d       = dictionary of p,q and seasonal components
/. r       > estimated parameters
i.estparam:{[endog;exog;errors;d;tr]
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
  if[tr;x:1f,'x];
  y:m#endog;
  first enlist[y]lsq flip x
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

// Predict a single ARMA/AR value
/* params   = fit parameters
/* dict     = dictionary containing indexing information
/* pvals    = list containing p lag values, q residual errors and previously predicted values
/* estres   = estimated residual error
/* typ      = typ of model bring used
/. r        > list of lag values, preds and predicted values
i.sngpred:{[params;exog;dict;pvals;estresid;typ]
  // entry in table from which prediction is being made
  exogval:exog count pvals 2;
  // check if trend is present
  pred:$[count[params]~count[m:exogval,get[".tm.i.",typ,"val"][pvals;dict]];
    // no trend value
    params mmu m;
    // add trend value to offset prediction (OoO not important as 2 single vectors)
    params[0]+m mmu 1_params];
  // if MA>0, estimate resid errors and append
  if[count pvals 1;estvals:exogval,$["SAR"~typ;#[neg[dict`n];];]pvals[0];
   pvals[1]:(1_pvals[1]),pred-mmu[estresid;estvals]];  
   // append new lag values, for next step calculations
  ((1_pvals[0]),pred;pvals[1];pvals[2],pred)}

// Extract fitted ARMA model params to return
/.r         > list of params needed for future predictions
i.ARMAparams:{[endog;coeff;dict;errs;n]
 (dict[`p]#neg[sum dict`q`p]#coeff;neg[dict`q]#coeff),
     (neg[n]#endog;neg[dict`q]#errs`err;errs`params),enlist dict
 }

// Extract fitted SARIMA model params to return
/.r         > list of params needed for future predictions
i.SARparams:{[endog;coeff;dict;errs;n]
 // number of seasonal components
 ns:count raze dict`P`Q;
 // Separate coeffs into normal and seasonal componants
 coefn:neg[ns]_coeff;coefs:neg[ns]#coeff;
 params:(dict[`p]#neg[sum dict`q`p]#coefn;neg[dict`q]#coefn;count[dict`P]#coefs;
   neg count[dict`Q]#coefs),((neg n|max[dict`P])#endog;
  (neg max raze dict`Q`q)#errs`err;errs`params);
 // Update dictionary values for seasonality funcs
 dict[`P`Q]:dict[`P`Q]-min each dict[`P`Q];
 params,enlist dict,enlist[`n]!enlist n
 }

// Extract appropriate lag and resid values for ARMA future predictions
/.r       > list of lag and resid values in appropriate order for predictions
i.ARMAval:{[pvals;dict] raze#[neg[dict`p];pvals[0]],pvals[1]}

// Extract appropriate lag and resid values for SARIMA future predictions
/.r      > list of lag and resid values in appropriate order for predictions
i.SARval:{[pvals;dict];
  raze#[neg[dict`p];pvals[0]],#[neg[dict`q];pvals[1]],
   pvals[0][dict[`P]],pvals[1][dict[`Q]]}

// Keys of ARMA and SARIMA model dictionaries
i.ARMAkeys:`params`tr_param`exog_param`p_param`q_param`lags`resid`estresid`pred_dict;
i.SARkeys:`params`tr_param`exog_param`p_param`q_param`P_param`Q_param`lags`resid`estresid`pred_dict;

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

