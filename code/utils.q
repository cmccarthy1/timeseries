// utility funcs for time series mdls

\d .tm

/*endog - variable to be predicted by model
/*exog - additional data table to be included when predicting data 
/*p - number of lags
/*d - number of differences to apply
/*q - number of residual errors
/*tr - include trend or not

// ARMA/AR models utils
// Fit ARMA mdl using hennan-Riessann 
/. r - model parameters used for predictions
i.ARMAmdl:{[endog;exog;p;q;tr]
 // convert exon table to matrix 
 if[98h~type exog;exog:"f"$i.mat[exog]];
 // build AR model to estimate resid errors
 estresid:ARfit[endog;exog;n:1+(p|q);0b];
 endogm:i.lagmat[endog;n];
 pred:((neg[count endogm]#exog),'endogm) mmu estresid`params;
 errors:(n _endog)-pred;
 // Using the resid errorrs calculate coefficients for ARMA model
 coef:i.estparam[endog;exog;errors;p;q;tr];
 // return dictionary with required values for forecasting
 `params`tr_param`exog_param`p_param`q_param`lags`resid`estresid!
   (coef;coef[tr-1];coef[tr+til count exog[0]];p#neg[q+p]#coef;
        neg[q]#coef;(neg n)#endog;(neg q)#errors;estresid`params)
 }

// Estimate parameters of ARMA mdl
/*errors - errors from AR mdl
/. r - estimated params for ARMA mdl
i.estparam:{[endog;exog;errors;p;q;tr]
 // create lag matrix for endow values
 endogm:i.lagmat[endog;p];
 // create lag matrix for resid errors
 resid:i.lagmat[errors;q];
 // decide how many values to use from training 
 m:neg count[endogm]&count[resid];
 // join exog, endow and resid values
 x:(m #exog),'(m #endogm),'m #resid;
 // add trend line if specified
 if[tr;x:1f,'x];
 // values to predict
 y:m #endog;
 // use least squared error method to get coefficient values
 inv[fx mmu x] mmu (fx:flip x) mmu y
 }

// Predict single ARMA/AR value
/* d - list of p lag values, q resid errors and the previous predicted values
/*estresid - The model params used for estimating resin errors
/ . r - list of p lag values, q resid errors and list of predicted values
i.sngpred:{[params;exog;p;d;estresid] 
 // check if trend is present
 pred:$[count[params]~count[m:exog[count d[2]],p _raze d[0 1]];
 // no trend, so multiply coefficients by values
 params mmu m;
 // add trend value
 params[0]+((1_params) mmu m)];
 // if MA>0, estimate error values
 errors:$[count d[1];pred-(estresid mmu exog[count d[2]],d[0]);()];
 // append new lag values, and resid errors for next step calculations
 ((1_d[0]),pred;(1_d[1]),errors;predn:d[2],pred)}


// Pure MA models utils

// Finding the params for pure MA models
/. r - returns innovations matrix to calculate past errors
i.innovations:{[data;q]
 v:(q+1)#0f;
 alpha:(q+1;q+1)#0f;
 v[0]:i.gamma[data;0];
 n:1;
  while[n<=q;
   k:0; 
   while[k<n;
   j:0; s:0f; 
    while[j<k;
    s+:(alpha[k;k-j]*alpha[n;n-j]*v[j]);
    j+:1; ];
   alpha[n;n-k]:(i.gamma[data;n-k] - s)%v[k]; 
   k+:1;
   s:0f; j:0;
   while[j<n;
   s+:(alpha[n;n-j] xexp 2)*v[j];
   j+:1; ];
  v[n]:v[0]-s];
  n+:1; ];
 alpha
 }

// Find residual errors for pure ma model
/*mat - innovations matrix
/. r - returns residual errors from model
i.residuals:{[data;q;mat]
 qdata:neg[q] #data;
 est:(q)#0f;
 n:1;
 while[n < q; 
        j:1;s:0f;
        while[j <= n;
            s:s+(mat[n;j]*(qdata[j-1]-est[j-1]));
            j+:1];
        est[n]:s;
        n+:1;
        s:0f
    ];
    qdata - est}

// AIC utils

// Fit the model with the params and score it using AIC
/*train - training data
/*test - testing data
/*param - list of parameters used (q,p,trend)
/. r - the AIC score to the corresponding params
i.fitscore:{[train;test;len;param]
 // create model using specified params
 mdl:ARIMAfit[train`endog;train`exog;param`p;param`d;param`q;param`tr];
 // predict values using model coeffs
 preds:ARIMApred[mdl;test`exog;len];
 // calculate aic score using predictions and true values
 i.aic[len# test`endog;preds;param]}

// Apply AIC scoring
/*true - true values
/*pred - predicted values
/. r - returns AIC score
i.aic:{[true;pred;params]
 // get residual sum of squares
 rss:sqr wsum sqr:true-pred;
 // get aic score
 sc:(2*k:sum params)+n*log(rss%n:count[pred]);
 //If k<40, use altered aic score
 $[k<40;sc+(2*k*(k-1))%n-k-1;sc]}

// Auto correlation function
/*h - order of autocorrelation
i.acf:{[data;h]
 i.gamma[h;data]%i.gamma[0;data]}


// General Utils

// Create a lagmatrix
/*data - data to create matrix 
/*lag - number of lags
i.lagmat:{[data;lag]
 n:count data;
 data til[n-lag]+\:til lag
 }

// Create matrix from table
i.mat:{flip value flip x}

/ Indicates if the data is stationary
/*r - returns boolean
i.stat:{[data]
 adfull:.ml.fresh.i.adfuller[data]`;
 $[adfull[1]<0.05;1b;0b]}

// Differencing to produce stationary time series
/d* - order of differencing
/ . r - list of differences in data
i.diff:{[data;d]
 d _d{deltas x}/data}

// AutoCorrelation function
/h* - order of autocorrelation
/. r -  autocorrelation of order h 
i.gamma:{[h;data]
 (neg[h]_data) cov h _data}


// Error calls

i.err.steps:{'`$"Exog length not long enough"} 
i.err.stat:{'`$"Time series not stationary, try another value of d"}
i.err.len:{'`$"Endog and Exog not the same length"}
