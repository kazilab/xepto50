import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from lmfit import Model, Parameters

# Logistic function from the Hill equation


def _4pl(log10_drug_con, slope, min_val, max_val, log10ic50):
    '''
    Args:
        log10_drug_con: Drug concentration in molar scale, log10 transformed.
        slope: Hill slope
        min_val: Minimum response
        max_val: Maximum response
        log10ic50: Log10 transformed IC50 in molar scale

    Returns:
        response: fitted response

        This is an implementation of the 4 Parameter Logistic (4PL) model, which is commonly used in bioassays and
        dose-response curve fitting. The 4PL model can describe sigmoidal curves, which are often observed in
        biological responses to varying doses of a drug or other substances.
        The parameters of the function:
        log10_drug_con: This parameter represents the logarithm base 10 of the drug concentration.
                        It is a variable on the x-axis in dose-response curves.
        slope: This parameter dictates the steepness of the curve in the sigmoidal 4PL model.
        min_val: This is the minimum value or the lower asymptote of the sigmoidal curve.
        max_val: This is the maximum value or the upper asymptote of the sigmoidal curve.
        log10ic50: This parameter is the logarithm base 10 of the IC50 value.
                   IC50 is a commonly used metric in pharmacology, representing the concentration of a substance
                   required to inhibit a biological process by half.
        The function calculates the response y, given the input parameters, using the 4PL equation:
            y = min_val + (max_val - min_val) / (1 + 10 ** (slope * (log10ic50 - log10_drug_con)))
        The result y is then returned by the function.
        This y value represents the response (on the y-axis of a dose-response curve) corresponding to the given drug
        concentration and the specified curve parameters. The 4PL model and the calculated y are useful for analyzing
        and visualizing the effect of a drug across different concentrations.
    '''

    y = min_val + (max_val - min_val) / (1 + 10 ** (slope * (log10ic50 - log10_drug_con)))
    return y


def fit_curve_4pl(log10con, inhibition):
    '''
    Args:
        log10con: Drug concentration in molar scale, log10 transformed.
        inhibition: A single column response values.

    Returns:
        min_inhibition: The minimum response.
        max_inhibition: The maximum rsponse.
        slope_curve_fit: Hill slope
        log10ic50_curve_fit: Log10 transformed IC50 in molar scale

        This function is designed to fit a set of data points to the 4 Parameter Logistic (4PL) model, making use of
        the _4pl function previously described. It takes as input:
        log10con: An array or list containing the logarithm base 10 of the drug concentrations.
        inhibition: An array or list containing the corresponding inhibition responses.
        The function performs the following operations:
            Initialize Parameters: It sets initial guesses for the parameters of the 4PL model, along with the lower
                and upper bounds for each parameter, to guide the curve fitting.
            Curve Fitting: The function employs the curve_fit method from scipy.optimize to estimate the parameters of
                the 4PL model that best fit the input data. The fitted parameters include the slope, minimum and
                maximum inhibitions, and log10ic50.
            Adjust Minimum and Maximum Inhibition Values: After obtaining the fitted parameters, it makes several
                adjustments to the minimum and maximum inhibition values to ensure they are within reasonable bounds
                (0-99 for minimum inhibition and 1-100 for maximum inhibition).
            Compute Running Mean: It computes the running average of the inhibition values and adjusts the maximum
                inhibition value based on this running average.
            Final Adjustments: It performs several checks and adjustments on the calculated parameters, specifically on
                the log10ic50 value, to ensure it falls within the range of the input log10 concentrations and makes
                sense given the inhibition values.
            Return Fitted Parameters: Finally, it returns the adjusted minimum inhibition, maximum inhibition, slope,
                and log10ic50 values that best describe the input data according to the 4PL model.

            This function is useful for bioassay and pharmacological studies where sigmoidal dose-response curves are
            analyzed to understand the effect of different drug concentrations on biological responses, represented by
            inhibition values in this case. By fitting the data to the 4PL model and obtaining key parameters,
            researchers can gain insights into the characteristics of the drug, such as its potency (IC50) and efficacy
            (maximal effect).
        Use curve_fit from scipy.optimize to estimate parameters, re-runs parameters and covariance
        Initial guesses: [slope, min_val, max_val, log10ic50]
    '''
    min_inhibition_ = np.min(np.nanmin(inhibition))
    max_inhibition_ = np.max(np.nanmax(inhibition))

    upper_min = max(min_inhibition_, 5)
    lower_max = min(max_inhibition_, 95)
    initial_guesses = [1, 0, 100, np.median(log10con)]
    lower_bounds = [0, 0, lower_max, np.min(log10con)]
    upper_bounds = [4, upper_min, 100, np.max(log10con)]

    curve_fit_results, _ = curve_fit(f=_4pl,
                                     xdata=log10con,
                                     ydata=inhibition,
                                     p0=initial_guesses,
                                     bounds=(lower_bounds, upper_bounds),
                                     maxfev=100000)
    # Extract estimated parameters by curve_fit
    slope_curve_fit, min_inhibition_curve_fit, max_inhibition_curve_fit, log10ic50_curve_fit = curve_fit_results
    # print('slope', slope_curve_fit, 'min_val', min_inhibition_curve_fit,
    #       'max_val', max_inhibition_curve_fit, 'logIC50', log10ic50_curve_fit)
    # put a control on min_val_ and min_val_est to set between 0 and 99
    min_inhibition_ori = max(0, min_inhibition_)
    min_inhibition_ori = min(99, min_inhibition_ori)

    min_inhibition_curve_fit = np.nanmax([0, min_inhibition_curve_fit])
    min_inhibition_curve_fit = np.nanmin([99, min_inhibition_curve_fit])
    min_inhibition = min(min_inhibition_ori, min_inhibition_curve_fit)
    # print('min_inhibition:', min_inhibition)

    # put a control on max_val and max_val_est to set between 0 and 100
    max_inhibition_ori = min(100, max_inhibition_)
    max_inhibition_ori = max(1, max_inhibition_ori)

    max_inhibition_curve_fit = np.nanmin([100, max_inhibition_curve_fit])
    max_inhibition_curve_fit = np.nanmax([1, max_inhibition_curve_fit])

    max_inhibition_pre = np.nanmax([max_inhibition_ori, max_inhibition_curve_fit])
    # print('max_inhibition_pre:', max_inhibition_pre)

    # Compute running mean
    running_average = pd.Series(inhibition).rolling(window=3).mean().to_numpy()
    # set the upper value of maximum value
    max_inhibition_run = max_inhibition_pre
    if np.nanmax(running_average[:-1]) > running_average[-1]:
        valid_mask = ~np.isnan(running_average)
        boolean_array = (running_average > running_average[-1]) & valid_mask
        max_inhibition_run = np.nanmax(inhibition[boolean_array])
    if np.nanmax(inhibition) > max_inhibition_run:
        valid_mask = ~np.isnan(inhibition)
        boolean_array = (inhibition > max_inhibition_run) & valid_mask
        max_inhibition_run = np.nanmean(inhibition[boolean_array]) + 1

    max_inhibition_post = np.nanmax([1, max_inhibition_run])
    max_inhibition_post = np.nanmin([100, max_inhibition_post])

    max_inhibition = np.nanmax([max_inhibition_pre, max_inhibition_post])
    # print(min_inhibition, max_inhibition, slope_curve_fit)

    # check if minimum and maximum values qualifies, otherwise set to the experimental values
    if min_inhibition >= max_inhibition:
        log10ic50_curve_fit = np.nanmax(log10con)
    elif log10ic50_curve_fit > np.nanmax(log10con):
        log10ic50_curve_fit = np.nanmax(log10con)
    elif log10ic50_curve_fit < np.nanmin(log10con):
        log10ic50_curve_fit = np.nanmin(log10con)

    if np.nanmax(inhibition) < 0:
        log10ic50_curve_fit = np.nanmax(log10con)
    elif np.nanmin(inhibition) > 100:
        log10ic50_curve_fit = np.nanmin(log10con)

    if np.nanmean(inhibition[:-2]) < 5:
        log10ic50_curve_fit = np.nanmax(log10con)

    return min_inhibition, max_inhibition, slope_curve_fit, log10ic50_curve_fit


def fit_lm_4pl(log10con, inhibition, curve_fitted_values):
    '''

    Args:
        log10con: Drug concentration in molar scale, log10 transformed.
        inhibition: A single column response values.
        curve_fitted_values: Outout from curve fitting.

    Returns:
        lmfit_fitted_model: lmfit model after fitting
        lm_fit_results: results from lmfit model

        This function utilizes the lmfit library to refine the curve fitting of the 4 Parameter Logistic (4PL) model
        on a given set of drug concentration (log10con) and corresponding inhibition values (inhibition). The function
        takes the initial estimates of the parameters (curve_fitted_values) as input, presumably from a previous
        fitting operation (e.g., the fit_curve_4pl function), and returns a more refined model as well as a tuple of
        relevant results.

        A step-by-step description of what this function does:
        Create lmfit Model: Initializes an lmfit Model object based on the previously defined _4pl function.
        Define Initial Parameters: Sets up the initial parameters for the model fitting based on the input
            curve_fitted_values and defined boundaries.
        First Fitting Attempt: Performs the first curve fitting attempt using the initial parameters and evaluates
            the model with the provided log10con and inhibition data.
        Second Fitting Attempt with Adjusted IC50: Adjusts the initial log10ic50 value to the median of log10con and
            performs a second curve fitting attempt.
        Decide on Fitting Result: Compares the standard deviation of residuals from both fitting attempts and chooses
            the model with lower residual standard deviation.
        Adjustment for Low Slope: If the slope value from the chosen model is less than or equal to 0.2, the function
            adjusts the slope bounds and the log10ic50 initial value, and performs another fitting.
        Extract Results and Errors: Extracts relevant results and error information from the fitted model.
        Return Fitted Model and Results: Returns the final fitted model and a tuple containing the minimum value,
            maximum value, error on log10ic50, and standard deviation of residuals.

        By following these steps, the fit_lm_4pl function refines the parameter estimates and achieves a better fit to
        the data using the lmfit library. The returned lmfit model object and the tuple of results can then be used for
        further analysis and interpretation of the biological responses to different drug concentrations.

    '''
    # Create a lmfit model
    min_inhibition, max_inhibition, initial_log10ic50, initial_slope = curve_fitted_values
    # print(min_inhibition, max_inhibition, initial_log10ic50, initial_slope)

    min_max_val = min(90, max_inhibition * 0.9)
    max_min_val = max(10, min_inhibition * 1.1)

    lmfit_model = Model(_4pl)

    # Define initial parameters
    lmfit_params = Parameters()
    lmfit_params.add('slope', value=initial_slope, min=0, max=4)
    # lmfit_params.add('max_val', value=max_val_est, min=lower_max, max=upper_max)
    lmfit_params.add('max_val', value=max_inhibition, min=min_max_val, max=100)
    lmfit_params.add('min_val', value=min_inhibition, min=0, max=max_min_val)
    lmfit_params.add('log10ic50', value=initial_log10ic50, min=np.min(log10con), max=np.max(log10con))

    # Fit the model
    lmfit_fitted_model1 = lmfit_model.fit(data=inhibition, params=lmfit_params, log10_drug_con=log10con)

    # Try with median as starting value for IC50
    lmfit_params['log10ic50'].value = np.median(log10con)
    lmfit_fitted_model2 = lmfit_model.fit(data=inhibition, params=lmfit_params, log10_drug_con=log10con)

    # Decide on which result to use
    ic50std_residual1 = np.std(lmfit_fitted_model1.residual)
    ic50std_residual2 = np.std(lmfit_fitted_model2.residual)
    lmfit_fitted_model = lmfit_fitted_model1 if ic50std_residual1 < ic50std_residual2 else lmfit_fitted_model2
    # print(result_ic50.params)
    # print(result_ic50.params['slope'].value)
    if lmfit_fitted_model.params['slope'].value <= 0.2:
        lmfit_params['log10ic50'].value = initial_log10ic50
        lmfit_params['slope'].min = 0.1
        lmfit_params['slope'].max = 2.5
        lmfit_fitted_model = lmfit_model.fit(data=inhibition, params=lmfit_params, log10_drug_con=log10con)

    error = lmfit_fitted_model.params['log10ic50'].stderr
    residual = np.std(lmfit_fitted_model.residual)
    min_val_lm = lmfit_fitted_model.params['min_val'].value
    max_val_lm = lmfit_fitted_model.params['max_val'].value
    lm_fit_results = min_val_lm, max_val_lm, error, residual
    return lmfit_fitted_model, lm_fit_results
