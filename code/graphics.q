\d .tm

// Load matplotlib
plt:.p.import[`matplotlib.pyplot]

// Plotting functionality parameters
/* data  = original dataset
/* vals  = calculated values
/* m     = bar indices
/* title = name of the graph

// Autocorrelation plot
acfplot:{[data]
  acf:i.acf[data;]each m:1_til 11&count[data];
  i.plotfn[data;acf;m;"AutoCorrelation"];}

// Partial Autocorrelation Plot
pacfplot:{[data]
  pacf:.ml.fresh.i.pacf[data;neg[1]+m:11&count data]`;
  i.plotfn[data;1_pacf;1_til m;"Partial AutoCorrelation"];}

// Plotting utility
i.plotfn:{[data;vals;m;title]
  conf:10#1.95%sqrt count data;
  plt[`:bar][m;vals;`width pykw 0.5];
  cfgkeys:`linewidth`linestyle`color`label;
  cfgvals:3,`dashed`red`conf_interval;
  plt[`:plot][m;conf;pykwargs cfgkeys!cfgvals];
  plt[`:legend][];
  plt[`:xlabel][`lags];
  plt[`:ylabel][`acf];
  plt[`:title][title];
  plt[`:show][];}
