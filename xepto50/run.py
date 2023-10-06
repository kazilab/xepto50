import pandas as pd
import numpy as np
import matplotlib.backends.backend_pdf
from .cal import XCal
from .utils import remove_row_outliers
import io


def x50(df,
        response='Viability',
        response_type='Percentage',
        drug_concentration_unit='Micromolar',
        remove_outlier=True,
        set_baseline=10,
        set_integration_limit=1,
        report_quality_scores=False
        ):
    '''

    Args:
        df: pandas dataframe
        response: Response upon drug, "Viability" or "Inhibition".
        response_type: How response was canculated "Percentage" or "Ratio".
        drug_concentration_unit: Concentration unit used.
        remove_outlier: Whether outlier to be removed.
        set_baseline: Set a baseline response.
        set_min_con: Set the minimum concentration.
        set_integration_limit: Set a intergation limit.

    Returns:
        final_df: Result df
        pdf_bytes: Figure

    '''
    # Get the list of existing columns
    existing_columns = df.columns.tolist()

    # Create a dictionary to rename the first four columns
    rename_dict = {
        existing_columns[0]: 'experiment_id',
        existing_columns[1]: 'cell_line',
        existing_columns[2]: 'drug_name',
        existing_columns[3]: 'drug_concentration'
    }

    # Rename the columns
    df_renamed = df.iloc[:, :4].rename(columns=rename_dict)

    # Get response columns in a DataFrame
    df_response = df.iloc[:, 4:]
    # print(df_response)
    # Check if single column or multiple columns
    num_columns = len(df_response.columns)
    if num_columns == 0:
        raise ValueError('Response data is missing')
    if response_type == 'Ratio' and np.nanmax(np.nanmax(df_response, axis=1)) <=5:
        df_response = df_response * 100
    elif np.nanmax(np.nanmax(df_response, axis=1)) <= 2:
        df_response = df_response * 100
    if response == 'Viability':
        df_response = 100 - df_response

    if num_columns == 1:
        response_mean = df_response
        response_sem = None
    elif num_columns == 2:
        response_mean = df_response.apply(np.nanmean, axis=1)
        response_sem = None
    else:
        if remove_outlier:
            max_iterations = 3
            previous_sem_sum = float('inf')
            df_imputed = pd.DataFrame()
            for i in range(max_iterations):

                # Remove outliers in each row
                df_outlier_remove = df_response.apply(remove_row_outliers, axis=1)

                # Use row mean to fill in missing values
                df_imputed = df_outlier_remove.apply(lambda row: row.fillna(np.nanmean(row)), axis=1)

                # Recalculate SEM after imputation
                response_sem_post = df_imputed.sem(axis=1)

                # Calculate sum of SEM
                current_sem_sum = response_sem_post.sum()

                if current_sem_sum < previous_sem_sum:
                    previous_sem_sum = current_sem_sum
                    df_response = df_imputed  # Set the imputed df as the new df
                else:
                    break  # Exit the loop if SEM is not reducing

            response_mean = df_imputed.apply(np.nanmean, axis=1)
            response_sem = df_imputed.sem(axis=1)
        else:
            response_mean = df_response.apply(np.nanmean, axis=1)
            response_sem = df_response.sem(axis=1)

    # print(response_sem)
    # df_response.to_excel('df_response.xlsx')
    df_renamed['inhibition_mean'] = response_mean
    if response_sem is not None:
        df_renamed['inhibition_sem'] = response_sem
    # print(df_renamed)
    experiment_ids = df_renamed['experiment_id'].unique()
    dfs = []
    pdf_buffer = io.BytesIO()
    pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_buffer)

    for exp_id in experiment_ids:
        df_id = df_renamed[df_renamed['experiment_id'] == exp_id]
        cell_lines = df_id['cell_line'].unique()
        for line in cell_lines:
            df_line = df_id[df_id['cell_line'] == line]
            drugs = df_line['drug_name'].unique()
            for drug in drugs:
                df_drug = df_line[df_line['drug_name'] == drug]
                if response_sem is not None:
                    drug_results, fig = XCal(drug_concentrations=df_drug['drug_concentration'],
                                             inhibition_values=df_drug['inhibition_mean'],
                                             inhibition_sem=df_drug['inhibition_sem'],
                                             concentration_unit=drug_concentration_unit,
                                             drug_name=drug,
                                             cell_line=line,
                                             experiment=exp_id,
                                             set_baseline_for_auc=set_baseline,
                                             set_integration_limit=set_integration_limit,
                                             report_quality_scores=report_quality_scores
                                             ).run()
                    dfs.append(drug_results)
                    pdf.savefig(fig)
                else:
                    drug_results, fig = XCal(drug_concentrations=df_drug['drug_concentration'],
                                             inhibition_values=df_drug['inhibition_mean'],
                                             inhibition_sem=None,
                                             concentration_unit=drug_concentration_unit,
                                             drug_name=drug,
                                             cell_line=line,
                                             experiment=exp_id,
                                             set_baseline_for_auc=set_baseline,
                                             set_integration_limit=set_integration_limit,
                                             report_quality_scores=report_quality_scores
                                             ).run()
                    dfs.append(drug_results)
                    pdf.savefig(fig)

    # Concatenate all dataframes
    final_df = pd.concat([pd.DataFrame([d]) for d in dfs], ignore_index=True)

    pdf.close()
    # Move the buffer position to the beginning of the stream
    pdf_buffer.seek(0)
    # Read the buffer content to bytes
    pdf_bytes = pdf_buffer.read()
    # Close the buffer
    pdf_buffer.close()
    return final_df, pdf_bytes


'''
# Save the PDF bytes to a file
with open('output_plots.pdf', 'wb') as f:
    f.write(pdf_bytes)
'''