import sys
import pandas as pd
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QPushButton, QFileDialog, QComboBox, QCheckBox,
                             QSpinBox, QLabel, QMainWindow, QAction, QGroupBox, QGridLayout, QSizePolicy, QHBoxLayout,
                             QDoubleSpinBox, QSizePolicy)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt, QThread, pyqtSignal

# Make sure to adjust the import statement according to your project structure
from .run import x50


class XeptoThread(QThread):
    result_signal = pyqtSignal('PyQt_PyObject')
    plots_signal = pyqtSignal('PyQt_PyObject')

    def __init__(self,
                 df,
                 response,
                 response_type,
                 drug_concentration_unit,
                 remove_outlier,
                 set_baseline,
                 set_integration_limit,
                 report_quality_scores
                 ):
        QThread.__init__(self)
        self.df = df
        self.response = response
        self.response_type = response_type
        self.drug_concentration_unit = drug_concentration_unit
        self.remove_outlier = remove_outlier
        self.set_baseline = set_baseline
        self.set_integration_limit = set_integration_limit
        self.report_quality_scores = report_quality_scores

    # def __del__(self):
    #    self.wait()

    def run(self):
        # simulate a long task
        result, plots = x50(df=self.df,
                            response=self.response,
                            response_type=self.response_type,
                            drug_concentration_unit=self.drug_concentration_unit,
                            remove_outlier=self.remove_outlier,
                            set_baseline=self.set_baseline,
                            set_integration_limit=self.set_integration_limit,
                            report_quality_scores=self.report_quality_scores
        )
        self.result_signal.emit(result)  # we use a signal to return the result
        self.plots_signal.emit(plots)


class Window(QMainWindow):
    def __init__(self):
        super(Window, self).__init__()

        # Set window title
        self.setWindowTitle("Xepto50")

        # Create menu bar
        menu = self.menuBar()

        # Create items for the menu bar
        file_menu = menu.addMenu("Menu")

        # Create actions for the file menu
        open_action = QAction("Load", self)
        open_action.triggered.connect(self.load_xlsx)  # connect to the load_csv method
        file_menu.addAction(open_action)

        process_action = QAction("Xepto", self)
        process_action.triggered.connect(self.process)  # connect to the process method
        file_menu.addAction(process_action)

        export_action = QAction("Export", self)
        export_action.triggered.connect(self.export_xlsx)  # connect to the export_csv method
        file_menu.addAction(export_action)

        export_action = QAction("Plots", self)
        export_action.triggered.connect(self.export_plots)  # connect to the export_csv method
        file_menu.addAction(export_action)

        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)  # connect to the close method
        file_menu.addAction(exit_action)

        # Create central widget
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        # Set fixed size
        #self.setFixedSize(600, 500)  # width, height

        # Grid layout
        self.layout = QGridLayout(self.central_widget)  # set layout on central widget


        # Set a title
        self.title = QLabel("Xepto50")
        self.title.setAlignment(Qt.AlignCenter)
        self.title.setStyleSheet("color: lightblue; font-size: 45px; font-weight: bold;")
        self.layout.addWidget(self.title, 0, 0, 1, 2)  # added to row 0, column 0

        self.kazilab = QLabel("KaziLab.se @ Lund University")
        self.kazilab.setAlignment(Qt.AlignRight)
        self.kazilab.setStyleSheet("color: #3399CC; font-size: 10px; font-weight: bold;")
        self.layout.addWidget(self.kazilab, 0, 2, 1, 2)  # added to row 0, column 0

        self.subtitle = QLabel("Calculate drug efficacy in batch mode!")
        self.subtitle.setAlignment(Qt.AlignRight)
        self.subtitle.setStyleSheet("color: #3399CC; font-size: 11px; font-weight: bold;")
        self.layout.addWidget(self.subtitle, 1, 2, 1, 2)  # added to row 0, column 0

        group_box0 = QGroupBox("Data to be used")  # the argument is the title of the group box
        font0 = QFont()
        font0.setPointSize(10)
        font0.setBold(True)
        group_box0.setFont(font0)
        group_box0.setStyleSheet("QGroupBox::title {color: #4a86e8;}")
        # Create a layout for the group box
        group_layout0 = QVBoxLayout()
        # Create an inner QGridLayout
        inner_layout0 = QGridLayout()


        # Label for load csv file button
        self.load_csv_label = QLabel("Load an Excel File:")
        self.load_csv_label.setAlignment(Qt.AlignRight)
        self.load_csv_label.setToolTip('Click "Load Excel file" button to load a data file.')
        self.load_csv_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        #self.layout.addWidget(self.load_csv_label, 1, 0, 1, 1)  # added to row 0, column 0

        # load csv file button
        self.load_csv_button = QPushButton("Load data")
        self.load_csv_button.setStyleSheet("QPushButton { font-size: 12px;}")
        self.load_csv_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.load_csv_button.setMinimumSize(120, 20)  # Set minimum size of QSpinBox (width, height)
        self.load_csv_button.setMaximumSize(120, 20)  # Set maximum size of QSpinBox (width, height)
        self.load_csv_button.setToolTip('Click "Load Excel file" button to load a data file.')
        self.load_csv_button.clicked.connect(self.load_xlsx)
        #self.layout.addWidget(self.load_csv_button, 1, 1, 1, 1)  # added to row 0, column 0, spans 1 row and 2 columns

        # label to display selected file path
        self.file_path_label = QLabel("")
        #self.layout.addWidget(self.file_path_label, 1, 2, 1, 2)  # added to row 1, column 0, spans 1 row and 2 columns

        # Shall we remove outliers
        self.remove_outlier_label = QLabel("Remove outliers:")
        self.remove_outlier_label.setAlignment(Qt.AlignRight)
        self.remove_outlier_label.setToolTip('Check this box if you want to remove outliers')
        self.remove_outlier_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        # self.layout.addWidget(self.mf_for_xgb_label, 3, 0, 1, 1)

        self.remove_outlier_box = QCheckBox()
        self.remove_outlier_box.setStyleSheet(
            "QCheckBox::indicator { width: 20px; height: 20px;} QCheckBox { font-size: 20px;}")
        # self.layout.addWidget(self.mf_for_xgb, 3, 1, 1, 1)
        self.remove_outlier_box.setToolTip('Check this box if you want to remove outliers')

        inner_layout0.addWidget(self.load_csv_label, 0, 0, 1, 1)
        inner_layout0.addWidget(self.load_csv_button, 0, 1, 1, 1)
        inner_layout0.addWidget(self.file_path_label, 0, 2, 1, 2)
        inner_layout0.addWidget(self.remove_outlier_label, 1, 0, 1, 1)
        inner_layout0.addWidget(self.remove_outlier_box, 1, 1, 1, 1)

        # Add inner_layout to the group_layout
        group_layout0.addLayout(inner_layout0)
        # Set the group box's layout
        group_box0.setLayout(group_layout0)
        # Add the group box to the main layout
        self.layout.addWidget(group_box0, 2, 0, 3, 4)

        group_box1 = QGroupBox("Set parameters")  # the argument is the title of the group box
        font1 = QFont()
        font1.setPointSize(10)
        font1.setBold(True)
        group_box1.setFont(font1)
        group_box1.setStyleSheet("QGroupBox::title {color: #4a86e8;}")
        # Create a layout for the group box
        group_layout1 = QVBoxLayout()
        # Create an inner QGridLayout
        inner_layout1 = QGridLayout()

        # response dropdown
        self.response_label = QLabel("Response type:")
        self.response_label.setAlignment(Qt.AlignRight)
        # self.response_label.setToolTip('Set a method to fill missing data with a preliminary number.')
        self.response_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        self.response_label.setToolTip('How drug response was measured, viability or inhibition?')
        # self.layout.addWidget(self.response_label, 2, 0, 1, 1)
        self.response_dropdown = QComboBox()
        self.response_dropdown.setStyleSheet("QComboBox { font-size: 12px;}")
        self.response_dropdown.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.response_dropdown.setMinimumSize(120, 20)  # Set minimum size of QSpinBox (width, height)
        self.response_dropdown.setMaximumSize(120, 20)  # Set maximum size of QSpinBox (width, height)
        self.response_dropdown.setToolTip('How drug response was measured, viability or inhibition?')
        self.response_dropdown.addItem("Viability")
        self.response_dropdown.addItem("Inhibition")
        # self.layout.addWidget(self.response_dropdown, 2, 1, 1, 1)  # added to row 1, column 1

        # response dropdown
        self.calculation_label = QLabel("Calculated as:")
        self.calculation_label.setAlignment(Qt.AlignRight)
        # self.calculation_label.setToolTip('Set a method to fill missing data with a preliminary number.')
        self.calculation_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        self.calculation_label.setToolTip('How drug response was calculated, as percentage or ratio?')
        # self.layout.addWidget(self.calculation_label, 2, 0, 1, 1)
        self.calculation_dropdown = QComboBox()
        self.calculation_dropdown.setStyleSheet("QComboBox { font-size: 12px;}")
        self.calculation_dropdown.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.calculation_dropdown.setMinimumSize(120, 20)  # Set minimum size of QSpinBox (width, height)
        self.calculation_dropdown.setMaximumSize(120, 20)  # Set maximum size of QSpinBox (width, height)
        self.calculation_dropdown.setToolTip('How drug response was calculated, as percentage or ratio?')
        self.calculation_dropdown.addItem("Percentage")
        self.calculation_dropdown.addItem("Ratio")
        # self.layout.addWidget(self.calculation_dropdown, 2, 1, 1, 1)  # added to row 1, column 1

        # response dropdown
        self.unit_label = QLabel("Concentration unit:")
        self.unit_label.setAlignment(Qt.AlignRight)
        # self.unit_label.setToolTip('Set a method to fill missing data with a preliminary number.')
        self.unit_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        self.unit_label.setToolTip('Unit of drug concentration?')
        # self.layout.addWidget(self.pre_imputation_label, 2, 0, 1, 1)
        self.unit_dropdown = QComboBox()
        self.unit_dropdown.setStyleSheet("QComboBox { font-size: 12px;}")
        self.unit_dropdown.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.unit_dropdown.setMinimumSize(120, 20)  # Set minimum size of QSpinBox (width, height)
        self.unit_dropdown.setMaximumSize(120, 20)  # Set maximum size of QSpinBox (width, height)
        self.unit_dropdown.setToolTip('Unit of drug concentration?')
        self.unit_dropdown.addItem("Nanomolar")
        self.unit_dropdown.addItem("Micromolar")
        self.unit_dropdown.addItem("Millimolar")
        self.unit_dropdown.addItem("Molar")
        self.unit_dropdown.addItem("Picomolar")
        # self.layout.addWidget(self.unit_dropdown, 2, 1, 1, 1)  # added to row 1, column 1
        # self.empty1_label = QLabel("")
        # self.empty2_label = QLabel("")
        # self.empty3_label = QLabel("")

        # set_integration_limit spinbox
        self.integration_label = QLabel("Integration limit:")
        self.integration_label.setAlignment(Qt.AlignRight)
        self.integration_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        self.integration_label.setToolTip(f'Set an integration limit to calculate AUC from log10IC50,'
                                          f'Value must be a log10 value, for example 0.1, 0.5, 1 etc.')

        self.integration_spinbox = QDoubleSpinBox()
        self.integration_spinbox.setStyleSheet("QDoubleSpinBox { font-size: 12px;}")
        self.integration_spinbox.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.integration_spinbox.setMinimumSize(110, 20)  # Set minimum size of QDoubleSpinBox (width, height)
        self.integration_spinbox.setMaximumSize(110, 20)  # Set maximum size of QDoubleSpinBox (width, height)
        self.integration_spinbox.setToolTip(f'Set an integration limit to calculate AUC from log10IC50,'
                                            f'Value must be a log10 value, for example 0.1, 0.5, 1 etc.')
        self.integration_spinbox.setRange(0.1, 2)
        self.integration_spinbox.setSingleStep(0.1)
        self.integration_spinbox.setDecimals(2)
        self.integration_spinbox.setValue(1)  # Setting default value to 1

        # Add some widgets to the group box's layout
        inner_layout1.addWidget(self.response_label, 0, 0, 1, 1)
        inner_layout1.addWidget(self.response_dropdown, 0, 1, 1, 1)
        # inner_layout.addWidget(self.empty1_label, 0, 2, 1, 2)

        inner_layout1.addWidget(self.calculation_label, 0, 2, 1, 1)
        inner_layout1.addWidget(self.calculation_dropdown, 0, 3, 1, 1)
        # inner_layout.addWidget(self.empty2_label, 1, 2, 1, 2)

        inner_layout1.addWidget(self.unit_label, 1, 0, 1, 1)
        inner_layout1.addWidget(self.unit_dropdown, 1, 1, 1, 1)
        inner_layout1.addWidget(self.integration_label, 1, 2, 1, 1)
        inner_layout1.addWidget(self.integration_spinbox, 1, 3, 1, 1)

        # Add inner_layout to the group_layout
        group_layout1.addLayout(inner_layout1)
        # Set the group box's layout
        group_box1.setLayout(group_layout1)
        # Add the group box to the main layout
        self.layout.addWidget(group_box1, 5, 0, 3, 4)

        group_box2 = QGroupBox("Set optional parameters")  # the argument is the title of the group box
        font2 = QFont()
        font2.setPointSize(10)
        font2.setBold(True)
        group_box2.setFont(font2)
        group_box2.setStyleSheet("QGroupBox::title {color: #4a86e8;}")
        # Create a layout for the group box
        group_layout2 = QVBoxLayout()
        # Create an inner QGridLayout
        inner_layout2 = QGridLayout()

        # set_baseline spinbox
        self.baseline_label = QLabel("Baseline response:")
        self.baseline_label.setAlignment(Qt.AlignRight)
        self.baseline_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        self.baseline_label.setToolTip(f'Set a baseline response (between 0 and 25, default value: 10)')

        self.baseline_spinbox = QSpinBox()
        self.baseline_spinbox.setStyleSheet("QDoubleSpinBox { font-size: 12px;}")
        self.baseline_spinbox.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.baseline_spinbox.setMinimumSize(110, 20)  # Set minimum size of QDoubleSpinBox (width, height)
        self.baseline_spinbox.setMaximumSize(110, 20)  # Set maximum size of QDoubleSpinBox (width, height)
        self.baseline_spinbox.setToolTip(f'Set a baseline response (between 0 and 25, default value: 10)')
        self.baseline_spinbox.setRange(0, 25)
        self.baseline_spinbox.setSingleStep(5)
        self.baseline_spinbox.setValue(0)  # Setting default value to 0


        # Shall we report quality scores
        self.report_quality_scores_label = QLabel("Fit quality scores:")
        self.report_quality_scores_label.setAlignment(Qt.AlignRight)
        self.report_quality_scores_label.setToolTip('Check this box if you want to report fir quality scores')
        self.report_quality_scores_label.setStyleSheet("QLabel {color: #3399CC; font-size: 12px; font-weight: bold;}")
        # self.layout.addWidget(self.mf_for_xgb_label, 3, 0, 1, 1)

        self.report_quality_scores_box = QCheckBox()
        self.report_quality_scores_box.setStyleSheet(
            "QCheckBox::indicator { width: 20px; height: 20px;} QCheckBox { font-size: 20px;}")
        # self.layout.addWidget(self.mf_for_xgb, 3, 1, 1, 1)
        self.report_quality_scores_box.setToolTip('Check this box if you want to report fir quality scores')

        self.status_label = QLabel("Status: Ready")
        self.status_label.setStyleSheet("color: #3399CC; font-size: 14px; font-weight: bold;")
        # self.layout.addWidget(self.status_label, 9, 2, 1, 2)  # You may want to adjust the position

        # Add some widgets to the group box's layout
        inner_layout2.addWidget(self.baseline_label, 0, 0, 1, 1)
        inner_layout2.addWidget(self.baseline_spinbox, 0, 1, 1, 1)
        # inner_layout.addWidget(self.empty1_label, 0, 2, 1, 2)

        inner_layout2.addWidget(self.report_quality_scores_label, 0, 2, 1, 1)
        inner_layout2.addWidget(self.report_quality_scores_box, 0, 3, 1, 1)

        inner_layout2.addWidget(self.status_label, 1, 1, 1, 2)

        # Add inner_layout to the group_layout
        group_layout2.addLayout(inner_layout2)
        # Set the group box's layout
        group_box2.setLayout(group_layout2)
        # Add the group box to the main layout
        self.layout.addWidget(group_box2, 8, 0, 3, 4)

        # button to handle the function
        self.process_button = QPushButton("Compute")
        self.process_button.clicked.connect(self.process)
        self.process_button.setToolTip('Press to run process - an input file must be provided to start the process')
        self.process_button.setFixedSize(200, 28)  # width, height
        self.process_button.setStyleSheet("QPushButton {font-size: 15px; font-weight: bold; color: #ffffff; background-color: #3399CC;}")
        container = QWidget()
        button_layout = QHBoxLayout()
        button_layout.addWidget(self.process_button, 0, Qt.AlignCenter)
        # Set QHBoxLayout to the container
        container.setLayout(button_layout)
        self.layout.addWidget(container, 12,0,1,2)


        # self.status_label = QLabel("Status: Ready")
        # self.status_label.setStyleSheet("color: #3399CC; font-size: 15px; font-weight: bold;")
        # self.layout.addWidget(self.status_label, 9, 2, 1, 2)  # You may want to adjust the position

        self.export_csv_button = QPushButton("Export result as XLSX")
        self.export_csv_button.clicked.connect(self.export_xlsx)
        self.export_csv_button.setToolTip('Press to export imputed data')
        self.export_csv_button.setFixedSize(200, 28)  # width, height
        self.export_csv_button.setStyleSheet("QPushButton {font-size: 15px; font-weight: bold; color: #ffffff; background-color: lightblue;}")
        container2 = QWidget()
        button_layout2 = QHBoxLayout()
        button_layout2.addWidget(self.export_csv_button, 0, Qt.AlignCenter)
        container2.setLayout(button_layout2)
        self.layout.addWidget(container2, 13,0,1,2)

        self.export_plots_button = QPushButton("Export plots as PDF")
        self.export_plots_button.clicked.connect(self.export_plots)
        self.export_plots_button.setToolTip('Press to export imputed data')
        self.export_plots_button.setFixedSize(200, 28)  # width, height
        self.export_plots_button.setStyleSheet("QPushButton {font-size: 15px; font-weight: bold; color: #ffffff; background-color: lightblue;}")
        container3 = QWidget()
        button_layout3 = QHBoxLayout()
        button_layout3.addWidget(self.export_plots_button, 0, Qt.AlignCenter)
        container3.setLayout(button_layout3)
        self.layout.addWidget(container3, 13,2,1,2)

        # Close button
        self.close_button = QPushButton("Close")
        self.close_button.clicked.connect(self.close)
        self.close_button.setFixedSize(200, 28)  # width, height
        self.close_button.setStyleSheet(
            "font-size: 15px; font-weight: bold; color: #ffffff; background-color: gray")
        container1 = QWidget()
        button_layout1 = QHBoxLayout()
        button_layout1.addWidget(self.close_button, 0, Qt.AlignCenter)
        container1.setLayout(button_layout1)
        self.layout.addWidget(container1, 12,2,1,2)  # You may want to adjust the position

    def load_xlsx(self):
        self.file_name, _ = QFileDialog.getOpenFileName(self, "Open XLSX", "", "XLSX Files (*.xlsx)")
        # update label text with selected file path
        self.file_path_label.setText(f"Loaded file: {self.file_name}")
        if self.file_name:
            print(f"Loaded file {self.file_name}")
        else:
            print("No file was selected")

    def process(self):
        if hasattr(self, 'file_name'):
            df = pd.read_excel(self.file_name)
            response = self.response_dropdown.currentText()
            response_type = self.calculation_dropdown.currentText()
            drug_concentration_unit = self.unit_dropdown.currentText()
            remove_outlier = self.remove_outlier_box.isChecked()
            set_baseline = self.baseline_spinbox.value()
            set_integration_limit = self.integration_spinbox.value()
            report_quality_scores = self.report_quality_scores_box.isChecked()

            self.status_label.setText("Running...")
            self.process_button.setEnabled(False)  # disable the process button while the task is running

            # start the task in a new thread
            self.thread = XeptoThread(df,
                                      response,
                                      response_type,
                                      drug_concentration_unit,
                                      remove_outlier,
                                      set_baseline,
                                      set_integration_limit,
                                      report_quality_scores
                                      )

            # connect signals
            self.thread.result_signal.connect(self.process_finished_result)
            self.thread.plots_signal.connect(self.process_finished_plots)
            self.thread.finished.connect(self.thread.deleteLater)
            self.thread.start()
            # self.thread.wait()

    def process_finished_result(self, result):
        self.status_label.setText("Completed - Ready to export results")
        self.process_button.setEnabled(True)  # re-enable the process button
        self.result = result  # save the result

    def process_finished_plots(self, plots):
        # self.status_label.setText("Completed - Ready to export plots")
        # self.process_button.setEnabled(True)  # re-enable the process button
        self.plots = plots  # save the result

    def export_xlsx(self):
        if hasattr(self, 'result'):
            # This will open a file dialog to select where to save the file
            path, _ = QFileDialog.getSaveFileName(self, "Save XLSX", "", "XLSX Files (*.xlsx)")
            if path:
                if not path.endswith('.xlsx'):
                    path += '.xlsx'
                # export the dataframe to a CSV file
                self.result.to_excel(path)
            else:
                print("No file location provided")
        else:
            print("No dataframe to export")

    def export_plots(self):
        if hasattr(self, 'result'):
            # This will open a file dialog to select where to save the file
            path, _ = QFileDialog.getSaveFileName(self, "Save plots", "", "PDF Files (*.pdf)")
            if path:
                if not path.endswith('.pdf'):
                    path += '.pdf'
                pdf_bytes = self.plots
                with open(path, 'wb') as f:
                    f.write(pdf_bytes)
            else:
                print("No file location provided")
        else:
            print("No plots to export")

if __name__ == "__main__":
    app = QApplication(sys.argv)

    window = Window()
    window.show()

    sys.exit(app.exec_())


def xgui():

    # Create a QApplication instance (needed for any GUI in PyQt)
    app = QApplication([])

    window = Window()
    window.show()

    # Start the application's event loop
    app.exec_()
