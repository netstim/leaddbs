from PyQt5.QtWidgets import QWidget
from functionalities import Functionalities
from pop_up_ui import pop_up_Mesh_Refinement


class PopUpMeshRefinement(Functionalities):
    def __init__(self, parent):
        self.main_win = QWidget()
        self.ui = pop_up_Mesh_Refinement.Ui_centralwidget()
        self.ui.setupUi(self.main_win)

        # Create main window object
        self.mainWindow = parent

        # set relative folder
        self.rel_folder = self.rel_folder()

        # Widget
        self.ui.widget_No_Of_Refinement_Frequencies.hide()
        self.ui.widget_No_Of_Refinement_Frequencies_2.hide()
        self.ui.widget_No_Of_Refinement_Frequencies_3.hide()

        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_No_Of_Refinement_Frequencies))
        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_No_Of_Refinement_Frequencies_2))
        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_No_Of_Refinement_Frequencies_3))

        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.ui.lineEdit_Refinement_Frequency.setText(
                "[-1]") if self.ui.comboBox_Refinement_Frequencies.currentText() == "No. of Refinement Frequencies" else None)

        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("No. of Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_Refinement_Frequencies))

        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("No. of Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_Refinement_Frequencies_2))

        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("No. of Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_Refinement_Frequencies_3))
        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.show_menu_item_on_combobox_state_change("No. of Refinement Frequencies",
                                                                 self.ui.comboBox_Refinement_Frequencies,
                                                                 self.ui.widget_Refinement_Frequencies_4))

        self.ui.comboBox_Refinement_Frequencies.currentIndexChanged.connect(
            lambda: self.ui.spinBox_No_Of_Refinement_Frequencies.setValue(
                -1) if self.ui.comboBox_Refinement_Frequencies.currentText() == "Refinement Frequencies" else None)

        # Save and cancel
        self.filename = "{}/pop_up_control/dictionaries/dict_mesh_refinement.py".format(self.rel_folder)

        self.ui.pushButton_Save.clicked.connect(
            lambda: self.saveCloseWindow(self.output_dict(), self.filename))
        self.ui.pushButton_Cancel.clicked.connect(lambda: self.closeWindow())

    def output_dict(self):
        output_dict = {
            'refinement_frequency': self.set_default_values([0],
                                                            self.ui.lineEdit_Refinement_Frequency.text()),
            'num_ref_freqs': self.set_default_values(-1,
                                                     self.ui.spinBox_No_Of_Refinement_Frequencies.value()),
            'rel_div_CSF': -1,
            'Adaptive_frac_div': self.set_default_values(0,
                                                         self.ui.doubleSpinBox_Deviation_Threshold.value()),
            'Min_Scaling': self.set_default_values(0, self.ui.doubleSpinBox_Min_Scaling.value()),
            'CSF_ref_reg': self.set_default_values(0, self.ui.doubleSpinBox_CSF_Ref_Reg.value()),
            'rel_div': self.set_default_values(0, self.ui.doubleSpinBox_Rel_Div.value()),
            'rel_div_current': self.set_default_values(0, self.ui.doubleSpinBox_Rel_Div_Current.value())
        }
        return output_dict

    def saveDict(self):
        self.saveCloseWindow(self.output_dict(), self.filename)
