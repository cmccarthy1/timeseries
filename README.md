# Example Time Series Notebooks

Time series is defined as a sequence of information collected sequentially in time. Unlike most supervised machine learning methods which takes every row in a data frame as an independent entry, time series forecasting takes into account trends observed in historical data in order to forecast future observations. This includes examining attributes such as seasonality components and lagged values among others. 

Throughout these notebooks various methods of time series forecasting are implemented on a range of datasets in order to compare implementation and results of models.

Both EmbedPy and JupyterQ are utilised in the forecasting of and plotting of results.

The following time series libraries are demonstrated in the notebooks:

1. [FBProphet](https://facebook.github.io/prophet/docs/quick_start.html): This is a time series model developed by Facebook for forecasting time series using an additive regression model where non-linear trends are fitted with seasonality trends such as daily, monthly yearly etc, along with holoday effects using the statistical model STAN.

2. [TBATS](https://github.com/intive-DataScience/tbats): This is a predictive model that uses exponential smoothing to forecast time serues with complex seasonal patterns.

3. [ARIMA](https://machinelearningmastery.com/arima-for-time-series-forecasting-with-python/): ARIMA is a auto regression model that uses past lagged values to forecast future values. Both a [statsmodel](https://www.statsmodels.org/stable/generated.statsmodel.tsa.arima_model.ARIMA.html) and kdb implementation are available.

4. [GluonTS](https://github.com/awslabs/gluon-ts): GluonTS utilises deep learning based methods in MXNet for probabilistic time series forecasting.

5. [LSTM](https://keras.io/layers/recurrent/): LSTM is an architrcture that utilises Recurrent Neural Networks (RNN) in order to learn contextual (lagged) information from sequential data.

In each notebook, two datasets are tested:

- Daily Temp: Daily Minimum temperatures reached in Melbourne over a 10 year period
- Bike Rental : The number of bikes rented every hour recorded by TFL

These two datasets were chosen to give a well rounded view of how to build both a simple and more complex model. 
`Daily Temp` represents a simple time series dataframe which requires very little data preparation and has no additional data columns to be included in the training of the model. `Bike Rental` however, is a more complex time series which requires additional data preparation and exogenous variables to be added to the models. 

**Notes**
These notebooks implement vanilla models of each time series forcasters with no feature extraction performed on the dataset. For a more in detailed description of these models and how to achieve best results, please follow the links provided in the notebooks.  

## Requirements

- kdb+>=? v3.5 64-bit
- Python 3.x
- [embedPy](https://github.com/KxSystems/embedPy)
- [JupyterQ](https://github.com/KxSystems/jupyterq)
- [ML-Toolkit](https://github.com/KxSystems/ml) (v0.3.x)

## Dependencies

Install the Python dependencies with

pip
```bash
pip install -r requirements.txt
```

or with conda
```bash
conda install --file requirements.txt
```
**N.B.** Additionally [graphviz](http://www.graphviz.org/download/) must be installed on the system running the notebooks.



