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

// Fit an ARMA model using the Hannan-Rissanen method
i.ARMAmdl:{[endog;exog;p;q;tr]
  // Convert exog to matrix
  if[98h=type exog;exog:"f"$i.mat[exog]];
  // Number of lags which must be accounted for
  n:1+p|q;
  // Construct an AR model to estimate the residual error parameters
  estresid:ARfit[endog;exog;n;0b]`params;
  // Convert the endogenous variable to lagged matrix
  endogm:i.lagmat[endog;n];
  // Predict future values based on estimations from AR model and use to estimate error
  err:(n _endog)-((neg[count endogm]#exog),'endogm)mmu estresid;
  // Based on the residual errors calculated the coefficients of the ARMA model
  coeff:i.estparam[endog;exog;err;p;q;tr];
  // Construct the dictionary for return from the function;
  key_vals:`params`tr_param`exog_param`p_param`q_param`lags`resid`estresid;
  params:(coeff(::;tr-1;tr+til count exog 0)),(p#neg[q+p]#coeff;neg[q]#coeff),
         (neg[n]#endog;neg[q]#err;estresid);
  key_vals!params
  }

// Estimate ARMA model parameters using OLS
i.estparam:{[endog;exog;errors;p;q;tr]
  // Create lagged matrices for the endogenous variable and residual errors
  endogm:i.lagmat[endog;p];
  resid :i.lagmat[errors;q];
  // Collect the data needed for estimation
  vals:(exog;endogm;resid);
  // How many data points are required
  m:neg min(count')1_vals;
  x:(,'/)m#'vals;
  // If required add a trend line variable
  if[tr;x:1f,'x];
  y:m#endog;
  first enlist[y]lsq flip x
  }

// Predict a single ARMA/AR value
/* params = fit parameters
/* exog   = exogenous variables
/* p      = number of parameters to ignore?
/* pvals  = list containing p lag values, q residual errors and previously predicted values
/* estres = estimated residual error
i.sngpred:{[params;exog;p;d;estresid]
  // entry in table from which prediction is being made
  exogval:exog count d 2;
  // check if trend is present
  pred:$[count[params]~count[m:exogval,p _raze d 0 1];
    // no trend value
    params mmu m;
    // add trend value to offset prediction (OoO not important as 2 single vectors)
    params[0]+m mmu 1_params];
  // if MA>0, estimate error values
  errors:$[count d 1;pred - estresid mmu exogval,d 0;()];
  // append new lag values, and resid errors for next step calculations
  ((1_d[0]),pred;(1_d[1]),errors;d[2],pred)}


// Akaike Information Criterion

/* true   = true values
/* pred   = predicted values
/* params = set of parameters used in creation of the model
i.aicscore:{[true;pred;params]
  // Calculate residual sum of squares, normalised for number of values
  rss:{wsum[x;x]%y}[true-pred;n:count pred];
  // Number of parameter
  k:sum params;
  aic:(2*k)+n*log rss;
  // if k<40 use the altered aic score
  $[k<40;aic+(2*k*k+1)%n-k-1;aic]
  }

// Fit a model, predict the test, return AIC score
/* train  = Training data as a dictionary with endog and exog data
/* test   = Testing data as a dictionary with endog and exog data
/* len    = Number of steps in the future to be predicted
/* params = parameters used in prediction
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
/. r > covariance between a dataset at time t and t-lag 
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


// Error calls
i.err.steps:{'`$"Exog length not long enough"}
i.err.stat:{'`$"Time series not stationary, try another value of d"}
i.err.len:{'`$"Endog and Exog not the same length"}

