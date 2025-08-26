# AnalysisSST
Analysis of sea surface temperature (SST) using time series models. Focusing on SST data from Year 1950 to Year 2022.

## Study objective
- Investigate the performance of the SGAR model on the analysis of SST behaviour in Labuan waters.
- Apply seasonal autoregressive integrated moving average (SARIMA) model on SST in Labuan waters.
- Apply neural network autoregressive (NNAR) model on SST in Labuan waters.
- Compare the forecast results provided by the models.

### Plots
- Increasing trend in SST
    <img src="images/time-series-plot-with-trend.png" width="200"/>
- Comparison with actual values
    | Model  | Plot |
    |--------|------|
    | SARIMA |![SARIMA](./images/sarima-forecast-vs-actual.png)|
    | NNAR   |![NNAR](./images/nnar-forecast-vs-actual.png)| 
    | SGAR   |![SGAR](./images/sgar-forecast-vs-actual.png)|
- Comparison of all results
    ![All results](./images/comparison-all-with-true.png)

## Acknowledgement
1. Dataset: Hadley Centre Sea Ice and Sea Surface Temperature data set ([HadISST](https://www.metoffice.gov.uk/hadobs/hadisst/)), in NetCDF format.
