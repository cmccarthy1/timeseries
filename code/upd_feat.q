// This is an updated version of the feat.q script providing some enhancements to the performance
// both of the lag and windowing functions, these functions are now less intertwinned and as such
// more friendly to integrations as standalone utilities to the ML toolkit

// Note the omission of any mention of forecasting points in these implementations, it is being
// left to the user to ensure their targets meet the criteria for not leaking data.

\d .ml

// Time series windowing function for equispaced time steps
/* tab       = simple table
/* col_names = symbol list of columns on which to apply windowing
/* funcs     = list of function names (as symbols) which are to be applied to the datasets
/* wins      = list of window sized on which to apply these functions
/. r         > table with functions applied on specified columns over appropriate windows
/.             remove the first max[wins] columns as these are produced with insufficient information
ts_window:{[tab;col_names;funcs;wins]
  // unique combinations of columns/windows and functions to be applied to the dataset
  uni_combs:(cross/)(funcs;wins;col_names);
  // column names for windowed functions (remove ".") to ensure that if namespaced columns
  // exist they don't jeopardize parsing of select statements.
  win_cols:`$ssr[;".";""]each sv["_"]each string uni_combs;
  // values from applied functions over associated windows
  win_vals:{i.swin[get string y 0;y 1;x y 2]}[tab]each uni_combs;
  max[wins]_tab,'flip win_cols!win_vals
  }

// Addition of lagged features to equi-spaced time-series data sets
/* tab       = simple table
/* col_names = names of the columns to lag
/* lags      = list of historic lags to be added as columns to the dataset 
/. r         > table with historical lags added, removing the first max[lags] rows 
/.             as these will be null due to insufficient historical data
ts_lag:{[tab;col_names;lags]
  if[1=count col_names;col_names,:()];
  if[1=count lags;lags,:()];
  lag_names:`$raze string[col_names],/:\:"_xprev_",/:string lags;
  lag_vals :raze xprev'[;tab col_names]each lags; 
  max[lags]_tab,'flip lag_names!lag_vals
  }


// Sliding window function
// Note: this is a modified version of that provided in qidioms, using floating point windows instead
//       of long windows to increase the diversity of functions that can be applied
/* f = function taking a single argument (vector) to be applied
/* w = window size to be applied
/* s = vector on which sliding window and associated function are to be applied
i.swin:{[f;w;s]f each{ 1_x,y }\[w#0f;s]}
