\d .tm

// Miscellaneous functionality relating to time-series analysis
// and model generation procedures


// @kind function
// @category misc
// @fileoverview Summary of the stationality of each vector of a multivariate time series 
//   or a single vector
// @param dset {dict/tab/num[]} a time series of interest, the entries should 
//   in each case be numeric data types.
// @return {keytab} informative outputs from the python adfuller test indicating
//   the stationality of each vector entry of the relevant dataset
stationality:{[dset]
  dtype:type dset;
  // Names to be provided to form the key for the return table
  keynames:$[99h=dtype;key dset;
    98h=dtype;cols dset;
    enlist`data
    ];
  // Column names associated with the returns from the augmented dickey fuller test
  dcols:`ADFstat`pvalue`stationary,`$raze each"CriticalValue_",/:string(1;5;10),\:"%";
  scores:i.stationaryScores[dset;dtype];
  keynames!flip dcols!scores
  }

// @kind function
// @category misc
// @fileoverview Retrieve the best parameters for an ARIMA model based on the
//   Akaike Information Criterion (AIC)
// @param train  {dict} training data dictionary containing `endog/`exog data
// @param test   {dict} testing data dictionary containing `endog/`exog data
// @param len    {integer} number of steps forward to predict
// @param params {dict} parameter sets to fit ARIMA model with 
// @return {dict} parameter set which produced the lowest AIC score
aicParam:{[train;test;len;params]
  i.dictCheck[;`endog`exog;]'[(train;test);("train";"test")];
  i.dictCheck[params;`p`d`q`tr;"params"];
  // get aic scores for each set of params
  scores:i.aicFitScore[train;test;len;]each flip params;
  // return best value
  bestScore:min scores;
  scoreEntry:enlist[`score]!enlist bestScore;
  params[;scores?bestScore],scoreEntry
  }


// Time-series feature engineering functionality

// @kind function
// @category misc
// @fileoverview Apply a set of user defined functions over variously sized sliding windows
//   to a subset of columns within a table
// @param tab      {tab} dataset onto which to apply the windowed functions
// @param colNames {symbol[]} names of the columns on which to apply the functions
// @param funcs    {symbol[]} names of the functions to be applied
// @param wins     {integer[]} list of sliding window sizes
// @return         {tab} table with functions applied on specified columns over 
//   appropriate windows remove the first max[wins] columns as these are produced
//   with insufficient information to be deemed accurate
tsWindow:{[tab;col_names;funcs;wins]
  // unique combinations of columns/windows and functions to be applied to the dataset
  uni_combs:(cross/)(funcs;wins;col_names);
  // column names for windowed functions (remove ".") to ensure that if namespaced columns
  // exist they don't jeopardize parsing of select statements.
  win_cols:`$ssr[;".";""]each sv["_"]each string uni_combs;
  // values from applied functions over associated windows
  win_vals:{i.slidingWindowFunction[get string y 0;y 1;x y 2]}[tab]each uni_combs;
  max[wins]_tab,'flip win_cols!win_vals
  }


// @kind function
// @category misc
// @fileoverview Apply a set of user defined functions over variously sized sliding windows
//   to a subset of columns within a table
// @param tab      {tab} dataset from which to generate lagged data
// @param colNames {symbol[]} names of the columns from which to retrieve lagged data
// @param lags     {integers[]} list of lagged values to retrieve from the dataset
// @return         {tab} table with columns added associated with the specied lagged
//   values, remove the first max[lags] rows as these are produced with insufficient 
//   information to be deemed accurate
tsLag:{[tab;col_names;lags]
  if[1=count col_names;col_names,:()];
  if[1=count lags;lags,:()];
  lag_names:`$raze string[col_names],/:\:"_xprev_",/:string lags;
  lag_vals :raze xprev'[;tab col_names]each lags;
  tab,'flip lag_names!lag_vals
  }


// Plotting functionality

// @kind function
// @category misc
// @fileoverview Plot and display an autocorrelation plot
// @param data {num[]} dataset from which to generate the autocorrelation plot
// @return {graph} display to standard out the autocorrelation bar plot
acfPlot:{[data]
  acf:i.autoCorrFunction[data;]each m:1_til 11&count[data];
  i.plotFunction[data;acf;m;"AutoCorrelation"];
  }

// @kind function
// @category misc
// @fileoverview Plot and display an autocorrelation plot
// @param data {num[]} dataset from which to generate the partial autocorrelation plot
// @return {graph} display to standard out the partial autocorrelation bar plot
pacfPlot:{[data]
  pacf:.ml.fresh.i.pacf[data;neg[1]+m:11&count data]`;
  i.plotFunction[data;1_pacf;1_til m;"Partial AutoCorrelation"];
  }

