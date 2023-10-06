import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
# Functions


def xplot(ic50,
          dose_interpolated,
          plot_min,
          plot_max,
          log10con,
          inhibition,
          inhibition_sem,
          log10dose_fine,
          fitted_curve_fine,
          conversion_unit,
          dose_interpolated_at,
          experiment,
          cell_line,
          drug_name,
          set_baseline_for_auc,
          set_integration_limit,
          ):
    """

    Args:
        ic50: Calculated IC50 value
        dose_interpolated: Dose interpolated, interpolated IC50 value
        plot_min: Min response value for plotting
        plot_max: Max response value for plotting
        log10con: Drug concentration log10 in molar
        inhibition: Response values
        inhibition_sem: SEM for response
        log10dose_fine: Log10 doses in molar for fitting, has 10000 values ranging from min log10con to max log10con
        fitted_curve_fine: Fitted response
        conversion_unit: Conversion unit
        dose_interpolated_at: Dose interpolated at concentration
        experiment: Experiment number
        cell_line: Cell line name
        drug_name: Drug name
        set_baseline_for_auc: Baseline response
        set_integration_limit: Set an integration limit for xepto50

    Returns:
        fig: matplotlib figure.

        This function, xplot, is designed to create a plot visualizing a dose-response curve, along with additional
        elements such as data points, interpolated values, IC50, areas under the curve (AUC), and annotations.
        The purpose of such a plot is to visualize how a drug’s inhibition changes with different concentrations and
        to show how well a fitted curve represents the actual data.

        Breakdown of the function:
            Parameters:
                ic50: Concentration of the drug that gives half-maximal response.
                dose_interpolated: The dose level at which interpolation is performed.
                plot_min, plot_max: Define the minimum and maximum for the y-axis scale.
                log10con: Log10 transformed drug concentrations.
                inhibition: Percentage inhibition values corresponding to log10con.
                inhibition_sem: Standard Error of the Mean for inhibition values.
                log10dose_fine: Fine-grained log10 transformed drug concentrations for plotting smooth curves.
                fitted_curve_fine: Fitted curve values corresponding to log10dose_fine.
                conversion_unit: Conversion factor for units of drug concentration.
                dose_interpolated_at: Percentage level at which interpolation is performed.
                experiment: Identifier for the experiment.
                cell_line: Name/Identifier of the cell line used.
                drug_name: Name of the drug.
                set_baseline_for_auc: Baseline value for AUC calculations.
                set_integration_limit: Integration limit for AUC calculations.
        Function Steps:
            Define Fitted Curve Function: Define a helper function fitted_curve that evaluates the lmfit_fitted_model
                at given dose ranges.
            Generate Data for Plotting: Create doses and curve_values arrays for plotting.
            Initialize Plot: Initialize a matplotlib figure and axis object.
            Set Y-Axis Limits and Plot Data Points: Set the limits for the y-axis and plot the actual data points,
                with error bars if inhibition_sem is provided.
            Plot Fitted Curves: Plot the fitted logistic curve and another fitted curve in different colors.
            Draw Vertical Lines and Shade AUC Areas: Draw vertical dashed lines at interpolated dose and log10IC50.
                Shade the area representing xAUC (extended AUC) in red and the regular AUC in light gray.
            Annotate Points: Annotate the interpolated dose and log10IC50 on the plot.
            Set Labels, Title, Legend, and Layout: Set the x-axis and y-axis labels, plot title, add a legend,
                adjust the layout, and close the plot.
            Return Figure: The function returns the matplotlib figure object containing the constructed plot.

        In summary, this function creates a comprehensive visual representation of a dose-response relationship,
        highlighting key parameters, fitted curves, and areas under the curves, which is useful for analyzing the
        efficacy and characteristics of a drug.

    """
    # Plotting
    min_scale = min(0, plot_min-2)
    max_scale = max(100, plot_max+2)
    log10dose_at_interpolated = np.log10(dose_interpolated * conversion_unit)
    log10ic50 = np.log10(ic50 * conversion_unit)

    fig, ax = plt.subplots()
    # Set y-axis limits
    ax.set_ylim([min_scale, max_scale])
    if inhibition_sem is not None:
        ax.errorbar(log10con, inhibition, yerr=inhibition_sem, fmt='o', color='red', label='Actual Data')
    else:
        ax.scatter(log10con, inhibition, color='red', label='Actual Data')
    ax.plot(log10dose_fine, fitted_curve_fine, color='lightblue', label='Fitted Logistic Curve')

    # Draw vertical lines at log10dose_at_50 and log10ic50_est
    ax.axvline(log10dose_at_interpolated, color='green', linestyle='--',
               label=f'Interpolated at {dose_interpolated_at}%')
    ax.axvline(log10ic50, color='blue', linestyle='--', label='IC50')
    """
    # Shade the area under the curve
    if np.max(np.max(fitted_curve_fine)) >=51:
        ax.fill_between(log10dose_fine,
                        fitted_curve_fine,
                        50,
                        where=((log10dose_fine >= log10dose_at_interpolated) &
                               (log10dose_fine <= log10dose_at_interpolated + set_integration_limit)),
                        color='lightgreen', label='xAUC')
    mask = fitted_curve_fine >= set_baseline_for_auc
    log10dose_fine_masked = log10dose_fine[mask]
    fitted_curve_fine_masked = fitted_curve_fine[mask]
    ax.fill_between(log10dose_fine_masked,
                    fitted_curve_fine_masked,
                    set_baseline_for_auc,
                    where=((log10dose_fine_masked >= np.min(log10con)) &
                           (log10dose_fine_masked <= np.max(log10con))),
                    color='lightgray', alpha=0.25, label='AUC')
    """
    # Annotate the points
    min_point = min_scale + 1
    mid_point = min_scale + 45
    ax.annotate(f'Interp@ {dose_interpolated_at}% = {dose_interpolated:.3f}',
                (log10dose_at_interpolated, mid_point), color='green', rotation=90)
    ax.annotate(f'IC50 = {ic50:.3f}', (log10ic50, min_point), color='purple', rotation=90)

    ax.set_xlabel('Log10 Concentration (Molar)')
    ax.set_ylabel('Inhibition %')
    ax.set_title(f'Curve for Exp: {experiment} Cell line: {cell_line} Drug: {drug_name}')
    ax.legend()
    plt.tight_layout()
    plt.close()
    return fig
