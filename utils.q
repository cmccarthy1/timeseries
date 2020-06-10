plt:.p.import[`matplotlib.pyplot]

// Include any missing datetimes in the table
/*dt - Datetime column name
/*tab -  table
/*tm - frequency of time in datetime col
datefill:{[dt;tab;tm]
         (flip enlist[dt]!enlist {x<max y}[;tab[dt]]{y+x}[tm]\min tab[dt])lj dt xkey tab}

// Train test split for time series
/*tab - input table
/*tar - target values
/* sz - train test split
ttstm:{[tab;tar;sz]`xtrain`ytrain`xtest`ytest!raze(tab;tar)@\:/:(0,floor n*1-sz)_til n:count tab}

// Plot the predicted and true values
/*preds - predicted values
/*true  -  true values
/*lab - plot label
pltresult:{
 [preds;true;lab]
 {plt[`:plot][x];}each(preds;true);
 plt[`:legend][`preds`true];
 plt[`:xlabel]["DateTime"];
 plt[`:ylabel][lab];
 plt[`:title][lab," vs DateTime"];
 plt[`:show][];}

// Plot the timeseries
/*dt - the timeframe of the timeseries
/*series - the values of the timeseries
/*lab - plot label
plttm:{[dt;series;lab]
 plt[`:plot][q2pydts dt;series];
 plt[`:xlabel]["Date"];
 plt[`:ylabel][lab];
 plt[`:title][lab," vs Date"];
 plt[`:show][];}

q2pydts:{.p.import[`numpy;
                   `:array;
                   "j"$x-("pmd"t)$1970.01m;
                   `dtype pykw "datetime64[",@[("ns";"M";"D");t:type[x]-12],"]"]}

// Reshape data to LSTM appropriate format
/*xdata - input data
/*ydata - target data
/*n_steps - number of time steps to use
/*n_feat - number of features being passed to the model
reshape:{[xdata;ydata;n_steps;n_feat]
  xdata:npa flip (m:count[xdata]-n)#'(-1_til n:n_steps+1)_\:xdata;
  xdata:xdata[`:reshape][m;n_steps;n_feat];
  `xdata`ydata!(xdata;(n)_ydata)}
