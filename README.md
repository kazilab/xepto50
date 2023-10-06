Xepto50 is developed for research purposes only.
The application has been tested using small datasets, please inform the writer if you find any problems.
The application has no guarantee and must be used with the user's responsibility.

We recommend using the Anaconda environment (install a recent copy of Anaconda from https://www.anaconda.com). Run Anaconda Powershell Prompt and use "pip install xputer" for installation.

To run the application GUI use "xepto50".

Jupyter Notebook script

from xepto50 import xgui

xgui()

data format: [Experiment] \t [Cell_line] \t [Drug_Name] \t [Drug_Concentration] \t [Response_Rep1] \t [Response_Rep2] \t [Response_Rep3] \t [Response_Rep4] \t [Response_Rep5] \t ........

Experiment: This column contains the experiment name or number.

Cell_line: This column indicates the name of the cell line or tissue.

Drug_Name: This column specifies the name of the drug.

Drug_concentration: This column shows the final drug concentration. The unit can be chosen from the GUI.

Response_Rep1: This column represents the drug response, which can be either viability or inhibition. It accepts values in percentage or ratio form. 

At least one response column is mandatory, but there is no specific limit on the number of response columns that can be added. If more than two response columns are provided, the system will calculate the SEM (Standard Error of the Mean). Users can also choose whether or not outliers should be removed.

-----------------------------
On Windows:
python setup.py install

On Linux:
sudo python setup.py install 

On Mac OS X:
sudo python ./setup.py install 

