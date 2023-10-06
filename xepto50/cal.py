import numpy as np
from .utils import *
import matplotlib
matplotlib.use('Agg')


class XCal:
    def __init__(self,
                 drug_concentrations: np.ndarray,
                 inhibition_values: np.ndarray,
                 inhibition_sem: np.ndarray = None,
                 concentration_unit: str = 'Nanomolar',
                 drug_name: str = None,
                 cell_line: str = None,
                 experiment: str = None,
                 set_baseline_for_auc: int = 10,
                 set_integration_limit: float = 1,
                 report_quality_scores: bool = False
                 )-> None:

        # initialize meta data
        self._initialize_meta_data(drug_name, cell_line, experiment, concentration_unit)
        # initialize data
        self._initialize_data(drug_concentrations, inhibition_values, inhibition_sem)
        # initialize parameters
        self._initialize_parameters(set_baseline_for_auc, set_integration_limit, report_quality_scores)
        # Updated later
        self._initialize_computational_properties()

    def _initialize_meta_data(self,
                              drug_name: str,
                              cell_line: str,
                              experiment: str,
                              concentration_unit: str
                              ) -> None:
        """
        Initializes metadata of the XCal object.

        :param drug_name: Name of the drug.
        :param cell_line: Name of the cell line.
        :param experiment: Name of the experiment.
        :param concentration_unit: Unit of the drug concentration.
        """
        self.drug_name = drug_name
        self.cell_line = cell_line
        self.experiment = experiment
        self.concentration_unit = concentration_unit

    def _initialize_data(self,
                         drug_concentrations: np.ndarray,
                         inhibition_values: np.ndarray,
                         inhibition_sem: np.ndarray
                         ) -> None:
        """
        Initializes data attributes of the XCal object.

        :param drug_concentrations: Array of drug concentrations.
        :param inhibition_values: Array of inhibition values.
        :param inhibition_sem: Array of standard error of mean for inhibition values.
        """
        # create a df with 'drug_concentration', 'inhibition', 'log10dose' and 'inhibition_sem
        self.df = create_df(drug_concentrations, inhibition_values, inhibition_sem, self.conversion_unit)
        self.drug_concentrations = self.df['drug_concentration'].values
        self.inhibition = self.df['inhibition'].values
        self.log10con = self.df['log10dose'].values
        self.min_log10dose = np.min(np.nanmin(self.log10con))
        self.max_log10dose = np.max(np.nanmax(self.log10con))
        if inhibition_sem is not None:
            self.inhibition_sem = self.df['inhibition_sem'].values
        else:
            self.inhibition_sem = None
        self.max_drug_concentration = np.nanmax(self.drug_concentrations)
        self.min_drug_concentration = np.nanmin(self.drug_concentrations)
        self.number_of_drug_concentrations = len(self.drug_concentrations)
        self.minimum_experimental_inhibition = np.nanmin(self.inhibition)
        self.maximum_experimental_inhibition = np.nanmax(self.inhibition)

    def _initialize_parameters(self,
                               set_baseline_for_auc: int,
                               set_integration_limit: float,
                               report_quality_scores: bool
                               ) -> None:
        """
        Initializes parameters of the XCal object.

        :param set_baseline_for_auc: Baseline value for AUC.
        :param set_integration_limit: Integration limit for calculations.
        :param report_quality_scores: Flag to report quality scores.
        """
        self.set_baseline_for_auc = set_baseline_for_auc
        self.set_integration_limit = set_integration_limit
        self.report_quality_scores = report_quality_scores

    def _initialize_computational_properties(self) -> None:
        """
        Initializes computational properties of the XCal object to None.
        """
        properties = [
            'lmfit_fitted_model', 'error', 'residual', 'log10dose_fine',
            'fitted_curve_fine', 'fitted_values', 'min_val_lm', 'max_val_lm',
            'min_inhibition', 'max_inhibition', 'dose_interpolated',
            'log10_dose_interpolated', 'dose_interpolated_at', 'ic50', 'auc',
            'xepto50', 'r2', 'ar2', 'syx', 'rmse', 'p', 'evar', 'merr', 'rmae',
            'mape', 'area', 'dss1', 'dss2', 'dss3', 'total_area'
        ]

        for prop in properties:
            setattr(self, prop, None)

    # convert concentration unit to value
    @property
    def conversion_unit(self) -> None:
        return unit_conversion(self.concentration_unit)

    def curve_fits(self):
        """
        Perform curve fitting and update relevant properties of the XCal instance.
        """
        self.initial_fit_results = self._initial_curve_fit()
        self._lm_fit(self.initial_fit_results)
        self._generate_fine_curves()
        self._adjust_min_max_inhibition()
        self._calculate_interpolated_values()
        self._calculate_ic50()
        pass

    def _initial_curve_fit(self):
        """
        Perform initial curve fitting using a 4PL model and return the results.
        """
        # use curve fit function to fit initial values
        min_inhibition, max_inhibition, log10ic50_curve_fit, slope_curve_fit = fit_curve_4pl(log10con=self.log10con,
                                                                                             inhibition=self.inhibition
                                                                                             )
        if min_inhibition == max_inhibition:
            max_inhibition += 0.001
        return min_inhibition, max_inhibition, log10ic50_curve_fit, slope_curve_fit


    def _lm_fit(self, initial_fit_results):
        """
        Improve the fitting using parameters from initial curve fitting.
        """
        self.lmfit_fitted_model, lm_fit_results = fit_lm_4pl(log10con=self.log10con,
                                                             inhibition=self.inhibition,
                                                             curve_fitted_values=initial_fit_results
                                                             )
        self.min_val_lm, self.max_val_lm, self.error, self.residual = lm_fit_results


    def _generate_fine_curves(self):
        """
        Generate a dense set of x values over the range of log10dose for a smoother curve.
        """
        self.log10dose_fine = np.linspace(self.min_log10dose, self.max_log10dose, 10000)
        self.fitted_curve_fine = self.lmfit_fitted_model.eval(log10_drug_con=self.log10dose_fine)
        self.fitted_values = self.lmfit_fitted_model.eval(log10_drug_con=self.log10con)
        # print(self.lmfit_fitted_model.params)

    def _adjust_min_max_inhibition(self):
        """
        Adjust minimum and maximum inhibition values.
        """
        min_inhibition = np.min([self.initial_fit_results[0], self.min_val_lm])
        self.min_inhibition = np.max([0, min_inhibition])
        max_inhibition = np.max([self.initial_fit_results[1], self.max_val_lm])
        self.max_inhibition = np.min([100, max_inhibition])


    def _calculate_interpolated_values(self):
        """
        Calculate interpolated values based on fitted curves.
        """
        # Find the index of the value in fitted_curve_fine closest to 50
        if self.min_inhibition <= 50 <= self.max_inhibition:
            index_at_interpolated = np.argmin(np.abs(self.fitted_curve_fine - 50))
            # Extract the corresponding log10dose value at this index
            self.dose_interpolated = ((10 ** self.log10dose_fine[index_at_interpolated]) /
                                      self.conversion_unit).round(2)
            self.dose_interpolated_at = self.fitted_curve_fine[index_at_interpolated].round(1)
        elif self.min_inhibition > 50:
            check_min_fitted_response = np.min(self.fitted_curve_fine)
            if check_min_fitted_response >= 50:
                self.dose_interpolated = ((10 ** np.min(self.log10dose_fine)) / self.conversion_unit).round(2)
                self.dose_interpolated_at = check_min_fitted_response.round(1)
            else:
                index_at_interpolated = np.argmin(np.abs(self.fitted_curve_fine - 50))
                # Extract the corresponding log10dose value at this index
                self.dose_interpolated = ((10 ** self.log10dose_fine[index_at_interpolated]) /
                                          self.conversion_unit).round(2)
                self.dose_interpolated_at = self.fitted_curve_fine[index_at_interpolated].round(1)
        elif self.max_inhibition < 50:
            check_max_fitted_response = np.max(self.fitted_curve_fine)
            if check_max_fitted_response < 50:
                self.dose_interpolated = ((10 ** np.max(self.log10dose_fine)) / self.conversion_unit).round(2)
                self.dose_interpolated_at = check_max_fitted_response.round(1)
            else:
                index_at_interpolated = np.argmin(np.abs(self.fitted_curve_fine - 50))
                self.dose_interpolated = ((10 ** self.log10dose_fine[index_at_interpolated]) /
                                          self.conversion_unit).round(2)
                self.dose_interpolated_at = self.fitted_curve_fine[index_at_interpolated].round(1)
        else:
            the_middle = (self.min_inhibition + self.max_inhibition) * 0.5
            # print(the_middle)
            index_at_interpolated = np.argmin(np.abs(self.fitted_curve_fine - the_middle))
            # Extract the corresponding log10dose value at this index
            self.dose_interpolated = ((10 ** self.log10dose_fine[index_at_interpolated]) /
                                      self.conversion_unit).round(2)
            self.dose_interpolated_at = self.fitted_curve_fine[index_at_interpolated].round(1)
        self.log10_dose_interpolated = np.log10(self.dose_interpolated * self.conversion_unit)


    def _calculate_ic50(self):
        """
        Calculate and adjust IC50 values.
        """
        # convert back to log10ic50 to IC50 of original unit
        self.ic50 = (10 ** self.lmfit_fitted_model.params['log10ic50'].value) / self.conversion_unit

        if self.lmfit_fitted_model.params['slope'].value < 0:
            self.ic50 = self.max_drug_concentration
        if self.lmfit_fitted_model.params['min_val'].value > 90:
            self.ic50 = self.min_drug_concentration

        if np.max([self.minimum_experimental_inhibition, self.maximum_experimental_inhibition]) < 10:
            self.ic50 = self.max_drug_concentration
        if self.minimum_experimental_inhibition > 50:
            self.ic50 = self.min_drug_concentration


    def scores(self):
        """
        Calculate and update various scores of the XCal instance.
        """
        self.auc, _ = cal_auc_aac(fitted_model=self.lmfit_fitted_model,
                                  set_baseline_for_auc=self.set_baseline_for_auc,
                                  max_inhibition=self.max_inhibition,
                                  log10con=self.log10con
                                  )
        self.xepto50 = xepto_score(lmfit_fitted_model=self.lmfit_fitted_model,
                                   dose_interpolated=self.log10_dose_interpolated,
                                   baseline_for_auc=self.set_baseline_for_auc,
                                   max_value_of_inhibition=self.max_inhibition,
                                   integration_limit=self.set_integration_limit
                                   )
        qs = quality_scores(original_values=self.inhibition,
                            fitted_values=self.fitted_values
                            )
        self.r2, self.ar2, self.syx, self.rmse, self.p, self.evar, self.merr, self.rmae, self.mape = qs

        for dss_type in range(1, 4):
            setattr(self, f'dss{dss_type}', self._calculate_dss(dss_type))


    def _calculate_dss(self, dss_type):
        """
        Calculate the DSS score based on the given type.
        """
        return dss(
            ic50=self.ic50,
            slope=self.lmfit_fitted_model.params['slope'].value,
            max_val=self.max_val_lm,
            min_con_tested=self.min_drug_concentration,
            max_con_tested=self.max_drug_concentration,
            y=self.set_baseline_for_auc,
            dss_type=dss_type,
            con_scale=self.conversion_unit
        )

    def generate_plot(self, plot_min, plot_max):
        """
        Generate the plot based on the instance attributes and the given min and max values.

        :param plot_min: Minimum value for plotting.
        :param plot_max: Maximum value for plotting.
        :return: A figure object.
        """
        return xplot(
            ic50=self.ic50,
            dose_interpolated=self.dose_interpolated,
            plot_min=plot_min,
            plot_max=plot_max,
            log10con=self.log10con,
            inhibition=self.inhibition,
            inhibition_sem=self.inhibition_sem,
            log10dose_fine=self.log10dose_fine,
            fitted_curve_fine=self.fitted_curve_fine,
            conversion_unit=self.conversion_unit,
            dose_interpolated_at=self.dose_interpolated_at,
            experiment=self.experiment,
            cell_line=self.cell_line,
            drug_name=self.drug_name,
            set_baseline_for_auc=self.set_baseline_for_auc,
            set_integration_limit=self.set_integration_limit
        )

    def run(self):
        """
        Execute the curve fitting and scoring methods, generate a plot, and compile the results.

        :return: A tuple containing the results dictionary and the figure object.
        """
        self.curve_fits()
        self.scores()
        plot_min = min(self.minimum_experimental_inhibition, self.min_val_lm)
        plot_max = max(self.maximum_experimental_inhibition, self.max_val_lm)
        fig = self.generate_plot(plot_min, plot_max)
        # Define common fields in drug_results
        drug_results = {
            'Experiment': self.experiment,
            'Cell line': self.cell_line,
            'Drug name': self.drug_name,
            'IC50': self.ic50.round(3),
            'IC interpolated': self.dose_interpolated,
            'IC interpolated at': self.dose_interpolated_at,
            'AUC': self.auc,
            'DSS1': self.dss1,
            'DSS2': self.dss2,
            'DSS3': self.dss3,
            'Xepto50': self.xepto50,
        }

        # Add additional fields based on the condition
        if self.report_quality_scores:
            additional_fields = {
                'Minimum drug concentration': self.min_drug_concentration,
                'Maximum drug concentration': self.max_drug_concentration,
                'Number of concentration': self.number_of_drug_concentrations,
                'Concentration unit': self.concentration_unit,
                'Minimum experimental inhibition': self.minimum_experimental_inhibition.round(2),
                'Minimum fitted inhibition': self.lmfit_fitted_model.params['min_val'].value.round(2),
                'Maximum experimental inhibition': self.maximum_experimental_inhibition.round(2),
                'Maximum fitted inhibition': self.max_val_lm.round(2) if self.max_val_lm != 100 else self.max_val_lm,
                'Slope': self.lmfit_fitted_model.params['slope'].value.round(2),
                'IC50 std error': self.error.round(3) if self.error is not None else None,
                'IC50 std residual': self.residual.round(3) if self.residual is not None else None,
                'R2': self.r2,
                'Adjusted R2': self.ar2,
                'RMSE': self.rmse,
                'Normality Test (p)': self.p,
                'Sy.x': self.syx,
                'Explained variance': self.evar,
                'Maximum residual error': self.merr,
                'Mean absolute error': self.rmae,
                'Mean abs percentage error': self.mape
            }
            drug_results.update(additional_fields)

        return drug_results, fig
