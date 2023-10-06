import numpy as np
from scipy.integrate import quad
from sklearn.metrics import (r2_score, mean_squared_error, explained_variance_score, max_error, mean_absolute_error,
                             mean_absolute_percentage_error)
from scipy import stats
import scipy.integrate as spi
from .curve_fitting import _4pl


def cal_auc_aac(fitted_model,
                set_baseline_for_auc,
                max_inhibition,
                log10con
                ):
    """
    Args:
        fitted_model: a model fitted for curve fitting, here we used lmfit
        set_baseline_for_auc: baseline response
        max_inhibition: maximum response
        log10con: drug concentrations in molar and log10 transformed

    Returns:
        auc: area under the curve
        aac: area above the curve

        The cal_auc_aac function calculates two metrics – Area Under the Curve (AUC) and Above the Curve (AAC) – from
        a fitted 4 Parameter Logistic (4PL) model. These metrics are often used in bioassay and pharmacological studies
        to quantify the overall behavior of the dose-response relationship.
        The function takes the following parameters as inputs:
            fitted_model: An lmfit Model object, presumably fitted using the fit_lm_4pl function described earlier.
            set_baseline_for_auc: The baseline value for calculating AUC.
            max_inhibition: The maximum inhibition value, one of the parameters of the 4PL model.
            log10con: An array or list containing the logarithm base 10 of the drug concentrations.

        The function proceeds with the following steps:
            Define Function for AUC: It defines a function fitted_function_auc that represents the fitted 4PL model
                adjusted by the provided baseline. This function will be used to calculate the AUC.
            Integrate for AUC: Using the quad function from scipy.integrate, it integrates fitted_function_auc over the
                range of log10con to compute the AUC – the area between the adjusted curve and the baseline.
            Define Function for AAC: Similarly, it defines another function fitted_function_aac to represent the maximum
                inhibition value minus the fitted 4PL model. This function will be used to calculate the AAC.
            Integrate for AAC: It integrates fitted_function_aac over the range of log10con to compute the
                AAC – the area above the curve up to the level of maximum inhibition.
            Round and Return Results: Finally, the function rounds the computed AUC and AAC to two decimal places and
            returns them as a tuple.

        By performing these steps, the cal_auc_aac function provides a quantitative assessment of the dose-response
        relationship depicted by the fitted 4PL model, which can be useful for comparing different drug effects or
        assessing drug efficacy and potency.
    """
    # Calculate AUC: Define the function using the fitted parameters, adjusted for the baseline

    def fitted_function_auc(dose_range, baseline):
        function = fitted_model.eval(log10_drug_con=dose_range) - baseline
        return function

    # Integrate over the range of interest to get the area between the curve and the baseline
    integration_result_auc = quad(fitted_function_auc, np.min(log10con), np.max(log10con),
                                  args=(set_baseline_for_auc,))

    # Calculate AAC: Define the function using the fitted parameters, adjusted for the baseline
    def fitted_function_aac(dose_2):
        function = max_inhibition - fitted_model.eval(log10_drug_con=dose_2)
        return function

    # Integrate over the range of interest to get the area between the curve and the baseline
    integration_result_aac = quad(fitted_function_aac, np.min(log10con), np.max(log10con))

    auc = round(integration_result_auc[0], 2)
    aac = round(integration_result_aac[0], 2)

    return auc, aac


def xepto_score(lmfit_fitted_model,
                dose_interpolated,
                baseline_for_auc,
                max_value_of_inhibition,
                integration_limit
                ):
    """
    Args:
        lmfit_fitted_model: Fitted lmfit model.
        dose_interpolated: Interpolated dose value
        baseline_for_auc: Baseline response
        max_value_of_inhibition: Maximum response
        integration_limit: Integration limit

    Returns:

        The xepto_score function calculates a specific score xepto50, based on a fitted model, presumably a 4 Parameter
        Logistic (4PL) model, from lmfit.
        This function is designed to quantify certain areas under the curve within specified limits and baseline
        adjustments.
        The inputs for the function are as follows:

            lmfit_fitted_model: The lmfit Model object that has been fitted with the data.
            dose_interpolated: The starting point for the integration.
            baseline_for_auc: The baseline value for calculating the area under the curve (AUC).
            max_value_of_inhibition: The upper limit used to calculate the area_of_interest.
            integration_limit: The range over which integration is to be performed.
        The function follows the subsequent steps to calculate the score:
            Define AUC Function: The function fitted_function_auc is defined using the parameters of the fitted model,
                adjusted by the baseline provided. It will be used to calculate the areas under the curve.
            Calculate Areas Under the Curve: The function integrates fitted_function_auc over the range specified
                by dose_interpolated and set_integration_limit. The integration_result_auc is adjusted by a fixed
                baseline of 50.
            Calculate Areas of Interest: The area_of_interest is calculated by multiplying set_integration_limit by the
                difference between max_value_of_inhibition and baseline_for_auc.
            Compute Xepto Score: The xepto50 score is computed by dividing the integration
            result by the area_of_interest and multiplying by 100, converting them into percentages.
            Return Results: The function rounds the computed xepto50.
        In essence, the xepto_score function offers a specialized assessment of the fitted model, producing score
        that quantify the proportion of the area under the curve within specific integration limits, with adjustments
        for baseline values. This score can be useful for comparing drug efficacy, potency, or the overall behavior
        of dose-response relationships in different scenarios.
    """
    # Calculate AUC: Define the function using the fitted parameters, adjusted for the baseline
    def fitted_function_auc(dose_range, baseline):
        function = lmfit_fitted_model.eval(log10_drug_con=dose_range) - baseline
        return function

    # Integrate over the range of log10IC50 to get the area between the curve and the baseline
    integration_result_auc = quad(fitted_function_auc, dose_interpolated,
                                    dose_interpolated + integration_limit, args=(50,))

    area_of_interest = (max_value_of_inhibition - baseline_for_auc) * integration_limit
    xepto50 = (integration_result_auc[0] / area_of_interest) * 100

    return round(xepto50, 2)


def quality_scores(original_values,
                   fitted_values
                   ):
    """
    Args:
        original_values: Original response
        fitted_values: Fitted response

    Returns:
        scores: quality control scores

        This quality_scores function is designed to calculate several statistical measures that quantify the quality of
        fit between a set of original (observed) values and a set of model-fitted (predicted) values.
        The function takes two parameters:
            original_values: A list or array containing the original observed values.
            fitted_values: A list or array containing the values predicted by a model.
        The function calculates the following quality scores and returns them as a tuple:
            R² Score (Coefficient of Determination): r2 measures the proportion of the variance in the dependent
                variable that is predictable from the independent variables. It ranges from 0 to 1, with 1 indicating
                perfect prediction.
            Adjusted R² Score: ar2 adjusts the R² score based on the number of observations and the number of
                predictors in the model. It is especially useful when comparing models with different numbers
                of predictors.
            Standard Error of the Estimate (Sy.x): syx provides a measure of the standard deviation of the errors
                (residuals) in the regression model.
            Root Mean Squared Error (RMSE): rmse is a frequently used measure of the differences between values
                predicted by a model and the values observed. It represents the square root of the second sample
                moment of the differences between predicted values and observed values or the quadratic mean of
                these differences.
            Shapiro-Wilk Normality Test P-value: The Shapiro-Wilk test tests the null hypothesis that the data was
                drawn from a normal distribution. A low p-value (< 0.05) indicates that the residuals have a
                distribution that is significantly different from a normal distribution.
            Explained Variance Score: evar measures the proportion to which a mathematical model accounts for the
                variation (dispersion) of a given data set.
            Maximum Residual Error: merr calculates the maximum difference between observed and predicted values.
            Root Mean Absolute Error (RMAE): rmae is the square root of the average of the absolute differences
                between the observed actual outcomes and the predictions made by the model.
            Mean Absolute Percentage Error (MAPE): mape expresses the forecast error as a percentage, which is useful
                when comparing errors across different scales.

        The function then returns a tuple containing all these calculated scores. Each of these scores provides a
        different perspective on the quality of the fit, helping to assess the performance and reliability of the
        regression model.
    """

    r2 = r2_score(fitted_values, original_values)
    # Adjusted R^2
    n = len(original_values)
    k = 1  # Number of predictors
    ar2 = 1 - (1 - r2) * (n - 1) / (n - k - 1)

    # Standard Error of the Estimate (Sy.x)
    syx = np.sqrt(sum((original_values - fitted_values) ** 2) / (n - 2))

    # RMSE
    rmse = np.sqrt(mean_squared_error(original_values, fitted_values))

    # Normality Test (Shapiro-Wilk)
    _, p = stats.shapiro(original_values - fitted_values)

    # Explained variance regression score function.
    evar = explained_variance_score(original_values, fitted_values)

    # The max_error metric calculates the maximum residual error.
    merr = max_error(original_values, fitted_values)

    # Mean absolute error regression loss.
    rmae = np.sqrt(mean_absolute_error(original_values, fitted_values))

    # Mean absolute percentage error (MAPE) regression loss.
    mape = mean_absolute_percentage_error(original_values, fitted_values)

    scores = r2, ar2, syx, rmse, p, evar, merr, rmae, mape

    return scores


def auc_cal(log10_drug_con_min,
            log10_drug_con_max,
            slope,
            set_baseline_for_auc,
            max_val,
            log10ic50
            ):
    """
    Args:
        log10_drug_con_min: Minimum drug concentration in log10 molar
        log10_drug_con_max: Maximum drug concentration in log10 molar
        slope: Hill slope
        set_baseline_for_auc: Baseline response
        max_val: Maximum response
        log10ic50: IC50 in log10 molar

    Returns:
        AUC: area under the curve

        This function auc_cal calculates the Area Under the Curve (AUC) for a 4 Parameter Logistic (4PL) model within a
        specified range of log10 drug concentrations. The 4PL model is often used to describe sigmoidal dose-response
        curves in pharmacology and bioassays. The AUC is a useful metric to quantify the overall effect of a drug.

        Here are the parameters of the function:
            log10_drug_con_min: The minimum log10 drug concentration for calculating the AUC.
            log10_drug_con_max: The maximum log10 drug concentration for calculating the AUC.
            slope: The slope of the curve at the inflection point.
            set_baseline_for_auc: The baseline (minimum response) value for the AUC calculation.
            max_val: The maximum response value of the curve.
            log10ic50: The log10 of the IC50 value, where IC50 is the concentration of the drug that gives half-maximal
                response.
        Within the function:
            The spi.quad function from Scipy is used to numerically integrate (calculate the area under) the 4PL curve
            between log10_drug_con_min and log10_drug_con_max.
            The _4pl function (which should be defined elsewhere in the code) calculates the value of the 4PL model at
            a given log10 drug concentration using the specified parameters (slope, set_baseline_for_auc, max_val,
            log10ic50).
            The quad function returns a tuple where the first element is the calculated AUC, and the second element is
            an estimate of the absolute error in the AUC. In this case, the function returns only the calculated AUC
            (result[0]).

            In summary, this function calculates and returns the AUC of a 4PL dose-response curve between specified
            log10 drug concentration limits using the provided curve parameters.
    """

    result = spi.quad(_4pl, log10_drug_con_min, log10_drug_con_max,
                      args=(slope, set_baseline_for_auc, max_val, log10ic50))

    return result[0]


def dss(ic50, slope, max_val, min_con_tested, max_con_tested, y=10, dss_type=2, con_scale=1e-9):
    a = float(max_val)
    b = float(slope)
    d = 0  # min response
    ic50 = float(ic50)
    min_con_tested = float(min_con_tested)
    max_con_tested = float(max_con_tested)
    min_con = np.log10(min_con_tested * con_scale)
    max_con = max_con_tested
    x2 = np.log10(max_con * con_scale)

    if (np.isnan(ic50) or np.isnan(b) or np.isnan(a) or
            np.isnan(min_con) or np.isnan(max_con)):
        return None

    if ic50 >= max_con:
        return 0

    if b == 0:
        return 0

    if a > 100:
        a = 100

    if b < 0:
        b = -b

    c = np.log10(ic50 * con_scale)

    # This is a logistic function used in Dotmatics.com
    # y = d+(a-d)/(1+10^(b*(c-x)))
    # inverse function
    # x = c - ((log(a-y)-log(d-y))/(b*log(10)))

    if a > y:
        if y != 0:
            x1 = c - ((np.log(a - y) - np.log(y - d)) / (b * np.log(10)))
            if x1 < min_con:
                x1 = min_con
            elif x1 > x2:
                x1 = x2
        else:
            x1 = min_con
        # the first part is integrate from x1 to x2 to get the area under the curve
        # int_y = ((((a - d) * np.log(1 + 10 ** (b * (c - x2)))) / (b * np.log(10)) + a * x2) -
        #         (((a - d) * np.log(1 + 10 ** (b * (c - x1)))) / (b * np.log(10)) + a * x1)) - y * (x2 - x1)

        area = spi.quad(_4pl, x1, x2, args=(slope, d, a, c))
        int_y = area[0]- y * (x2 - x1)

        total_area = (x2 - min_con) * (100 - y)
        norm_area = 0
        if dss_type == 1:
            norm_area = (int_y / total_area) * 100
        elif dss_type == 2:
            norm_area = (int_y / total_area) * 100 / np.log10(a)
            if norm_area > 50:
                norm_area = 0
        elif dss_type == 3:
            norm_area = ((int_y / total_area) * 100 * (np.log10(100) / np.log10(a)) *
                         ((x2 - x1) / (x2 - min_con)))

        if norm_area < 0 or norm_area > 100:
            return 0
        else:
            return round(norm_area, 2)

    else:
        return 0