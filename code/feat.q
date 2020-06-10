// Feature Engineering time series functions

\d .tm

/*tab - table to create features from
/*fp - forecast point( the time a prediction is being made)
/*col - list of columns to apply the transformations to

// Create window feature
/*n - list of window sizes over which to apply the functions
/*fnc - list of symbol functions to apply over the window (eg`min `max `avg)
/r. - returns the table with the additional window feature columns
addFeatWind:{[tab;n;fp;col;fnc]i.addFeat[tab;n;fp;col;fnc;"wnd"]}

// Create Lagged value features
/*n -  list of past values (from forecast point) to extract
/r. - returns the table with the additional lagged feature columns
addFeatLag:{[tab;n;fp;col]i.addFeat[tab;n;fp;col;`xprev;"lag"]}

// Update table with new feature columns
/*typ - the type of feature extraction (lagged (`prev) / windowed (`wnd))
i.addFeat:{[tab;n;fp;col;fnc;typ]![tab;();0b;raze i.applyFunc[fnc cross n;;typ;fp]each col]}

// Create dictionary of new columns and corresponding functions to extract the data
/*coln - single name of column to apply feature to
/r. - returns a dictionary 
i.applyFunc:{[fnc;coln;typ;fp](`$"_"sv/:string coln,'fnc)!
        {[x;y;z;fp]get[".tm.i.fnc",z][y;x[1];get string x[0];fp]}[;coln;typ;fp]each fnc}

// Windowed features being applied to a column
i.fncwnd:{[x;y;z;fp](,;fp#0N;(_;neg[fp];(`.tm.i.windowfnc;x;y;z)))}

// Lagged features being applied to a column
i.fnclag:{[x;y;z;fp](,;fp#0N;(_;neg[fp];(z;y;x)))}

// Create windowed in dataset and apply function
i.windowfnc:{z each x[neg[y-1]+til[y]+/:til count x]}
