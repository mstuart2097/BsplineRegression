# Price Yield Correlation and Revenue Insurance Premiums
The code in this repository provides details for performing our semi-parametric quantile regression (SQR) with penalized B-splines to assess if the correlation between price and yield and the revenue insurance premiums for both corn and soybeans from the years 1990-2023 depends on the level of leftover stocks/storage, defined as leftover yield from previous harvests.  

The results can be found in our manuscript titled "The Impact of Stocks on Correlations between Crop Yields and Prices and on Revenue Insurance Premiums using Semiparametric Quantile Regression".

The code files are labelled as follows:
- BSplineFunctions are the code to calculate the appropriate B-Spline traditional regression and quantile regression coefficients for a given quantile level and smoothing parameter for the penalty term
- Sim... are the code for conducting the simulation study as well as creating the plots for our simulation study
- Corn_by_State_1990_pbyt is the code for running the SQR by time and state for each quantile regression (Note: The associated file starting with Soybeans is the same code for the soybean dataset)
- Corn_Plots_State_1990_pbyt is the code for creating the plots from the SQR output (Note: The associated file starting with Soybeans is the same code for the soybean dataset)

The data files are found directly from the United States Department of Agriculture's National Agricultural Statistics Service database at https://quickstats.nass.usda.gov/.  The information on GDP is found using end-of-year GDP data from the Federal Reserve database.  The information on implied volatility was found from February futures prices for November's harvest for Corn and March futures prices for November's harvest for Soybeans.
