\d .tm

plt:.p.import[`matplotlib.pyplot]

// Autocorrelation plot
acfplot:{[data]
 acfvals:.tm.i.acf[data;]each m:1_til 11&count[data];
 conf:10#1.95%sqrt count data;
 plt[`:bar][m;acfvals;`width pykw 0.5];
 plt[`:plot][m;conf;`linestyle pykw `dashed;`color pykw "red";`linewidth pykw 3;`label pykw "conf_interval"];
 plt[`:legend][];
 plt[`:xlabel][`lags];
 plt[`:ylabel][`acf];
 plt[`:title][`Autocorrelation];
 plt[`:show][];}

// Partial Autocorrelation Plot
pacfplot:{[data]
 pacfvals:.ml.fresh.i.pacf[data;neg[1]+m:11&count data]`;
 conf:10#1.95%sqrt count data;
 plt[`:bar][1_til m;1_pacfvals;`width pykw 0.5];
 plt[`:plot][1_til m;conf;`linestyle pykw `dashed;`color pykw "red";`linewidth pykw 3;`label pykw "conf_interval"];
 plt[`:legend][];
 plt[`:xlabel][`lags];
 plt[`:ylabel][`acf];
 plt[`:title]["Partial AutoCorrelation"];
 plt[`:show][]; }
