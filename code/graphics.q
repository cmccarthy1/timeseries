\d .tm

plt:.p.import[`matplotlib]`:pyplot

// Autocorrelation plot
acfplot:{[data]
 acfvals:i.acf[data;]each m:1_til 11&count[data];
 plt[`:bar][m;acfvals];
 plt[`:xlabel][`lags];
 plt[`:ylabel][`acf];
 plt[`:title][`Autocorrelation];
 plt[`:show][];}

// Partial Autocorrelation Plot
pacfplot:{[data]
 pacfvals:.ml.fresh.i.pacf[data;neg[1]+m:11&count data]`;
 plt[`:bar][1_til m;1_pacfvals];
 plt[`:xlabel][`lags];
 plt[`:ylabel][`acf];
 plt[`:title]["Partial AutoCorrelation"];
 plt[`:show][]; }
