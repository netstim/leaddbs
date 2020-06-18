from PyQt5.QtWidgets import QWidget
from functionalities import Functionalities
from pop_up_ui import pop_up_External_Neuron_Array


class PopUpExternalNeuronArray(Functionalities):
    def __init__(self, parent):
        self.main_win = QWidget()
        self.ui = pop_up_External_Neuron_Array.Ui_centralwidget()
        self.ui.setupUi(self.main_win)

        # Create main window object
        self.mainWindow = parent

        # set relative folder
        self.rel_folder = self.rel_folder()

        # initial states
        self.ui.checkBox_Global_Rot.setCheckState(2)
        self.ui.widget_Global_Rot_Unchecked.hide()

        self.ui.checkBox_Global_Rot.stateChanged.connect(
            lambda: self.enable_disable_menu_item_on_checkbok_click_double_menu(self.ui.checkBox_Global_Rot,
                                                                                self.ui.widget_Global_Rot_Checked,
                                                                                self.ui.widget_Global_Rot_Unchecked))

        self.ui.lineEdit_X_Coordinate_Old.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_X_Coordinate_Old.text(),
                                                        self.ui.lineEdit_X_Coordinate_Old))
        self.ui.lineEdit_Y_Coordinate_Old.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_Y_Coordinate_Old.text(),
                                                        self.ui.lineEdit_Y_Coordinate_Old))
        self.ui.lineEdit_Z_Coordinate_Old.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_Z_Coordinate_Old.text(),
                                                        self.ui.lineEdit_Z_Coordinate_Old))
        self.ui.lineEdit_XY_Angles.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_XY_Angles.text(), self.ui.lineEdit_XY_Angles))
        self.ui.lineEdit_YZ_Angles.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_YZ_Angles.text(), self.ui.lineEdit_YZ_Angles))
        self.ui.lineEdit_ZX_Angles.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_ZX_Angles.text(), self.ui.lineEdit_ZX_Angles))
        self.ui.lineEdit_Alpha_Array_Glob.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_Alpha_Array_Glob.text(),
                                                        self.ui.lineEdit_Alpha_Array_Glob))
        self.ui.lineEdit_Beta_Array_Glob.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_Beta_Array_Glob.text(),
                                                        self.ui.lineEdit_Beta_Array_Glob))
        self.ui.lineEdit_Gamma_Array_Glob.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered(self.ui.lineEdit_Gamma_Array_Glob.text(),
                                                        self.ui.lineEdit_Gamma_Array_Glob))

        # Save and cancel
        self.filename = "{}/pop_up_control/dictionaries/dict_external_neuron_array.py".format(self.rel_folder)

        self.ui.pushButton_Save.clicked.connect(
            lambda: self.saveCloseWindow(self.output_dict(), self.filename))
        self.ui.pushButton_Cancel.clicked.connect(lambda: self.closeWindow())

    def output_dict(self):
        if self.ui.checkBox_Global_Rot.checkState() == 2:
            output_dict = {
                'Global_rot': self.corrector(1, self.ui.checkBox_Global_Rot.checkState()),
                'x_seed': self.set_default_values(0, self.ui.doubleSpinBox_X_Seed.value()),
                'y_seed': self.set_default_values(0, self.ui.doubleSpinBox_Y_Seed.value()),
                'z_seed': self.set_default_values(0, self.ui.doubleSpinBox_Z_Seed.value()),
                'x_steps': self.set_default_values(0, self.ui.spinBox_X_Steps.value()),
                'y_steps': self.set_default_values(0, self.ui.spinBox_Y_Steps.value()),
                'z_steps': self.set_default_values(0, self.ui.spinBox_Z_Steps.value()),
                'x_step': self.set_default_values(0, self.ui.doubleSpinBox_X_Step.value()),
                'y_step': self.set_default_values(0, self.ui.doubleSpinBox_Y_Step.value()),
                'z_step': self.set_default_values(0, self.ui.doubleSpinBox_Z_Step.value()),
                'alpha_array_glob': self.set_default_values([0],
                                                            self.ui.lineEdit_Alpha_Array_Glob.text()),
                'beta_array_glob': self.set_default_values([0],
                                                           self.ui.lineEdit_Beta_Array_Glob.text()),
                'gamma_array_glob': self.set_default_values([0],
                                                            self.ui.lineEdit_Gamma_Array_Glob.text())
            }
        else:
            output_dict = {
                'Global_rot': self.corrector(1, self.ui.checkBox_Global_Rot.checkState()),
                'X_coord_old': self.set_default_values(0, self.ui.lineEdit_X_Coordinate_Old.text()),
                'Y_coord_old': self.set_default_values(0, self.ui.lineEdit_Y_Coordinate_Old.text()),
                'Z_coord_old': self.set_default_values(0, self.ui.lineEdit_Z_Coordinate_Old.text()),
                'YZ_angles': self.set_default_values([0], self.ui.lineEdit_YZ_Angles.text()),
                'ZX_angles': self.set_default_values([0], self.ui.lineEdit_ZX_Angles.text()),
                'XY_angles': self.set_default_values([0], self.ui.lineEdit_XY_Angles.text())
            }

        return output_dict

    def saveDict(self):
        self.saveCloseWindow(self.output_dict(), self.filename)
