from PyQt5.QtWidgets import QWidget
from functionalities import Functionalities
from pop_up_ui import pop_up_Full_Field_IFFT


class PopUpFullFieldIFFT(Functionalities):
    def __init__(self, parent):
        self.main_win = QWidget()
        self.ui = pop_up_Full_Field_IFFT.Ui_centralwidget()
        self.ui.setupUi(self.main_win)

        # Create main window object
        self.mainWindow = parent

        # set relative folder
        self.rel_folder = self.rel_folder()


        self.ui.checkBox_VTA_From_DivE.stateChanged.connect(
            lambda: self.connected_check_boxes_3(self.ui.checkBox_VTA_From_DivE, self.ui.checkBox_VTA_From_NEURON,
                                                 self.ui.checkBox_VTA_From_E))
        self.ui.checkBox_VTA_From_NEURON.stateChanged.connect(
            lambda: self.connected_check_boxes_3(self.ui.checkBox_VTA_From_NEURON, self.ui.checkBox_VTA_From_DivE,
                                                 self.ui.checkBox_VTA_From_E))
        self.ui.checkBox_VTA_From_E.stateChanged.connect(
            lambda: self.connected_check_boxes_3(self.ui.checkBox_VTA_From_E, self.ui.checkBox_VTA_From_NEURON,
                                                 self.ui.checkBox_VTA_From_DivE))

        self.ui.checkBox_VTA_From_DivE.stateChanged.connect(
            lambda: self.show_menu_item_on_checkbox_click(self.ui.checkBox_VTA_From_DivE,
                                                          self.ui.widget_Activation_Threshold_1))

        self.ui.checkBox_VTA_From_DivE.stateChanged.connect(
            lambda: self.show_menu_item_on_checkbox_click(self.ui.checkBox_VTA_From_DivE,
                                                          self.ui.widget_Activation_Threshold_2))

        self.ui.checkBox_VTA_From_DivE.stateChanged.connect(
            lambda: self.show_menu_item_on_checkbox_click(self.ui.checkBox_VTA_From_DivE,
                                                          self.ui.widget_Activation_Threshold_3))

        self.ui.checkBox_VTA_From_E.stateChanged.connect(
            lambda: self.show_menu_item_on_checkbox_click(self.ui.checkBox_VTA_From_E,
                                                          self.ui.widget_Activation_Threshold_1))
        self.ui.checkBox_VTA_From_E.stateChanged.connect(
            lambda: self.show_menu_item_on_checkbox_click(self.ui.checkBox_VTA_From_E,
                                                          self.ui.widget_Activation_Threshold_2))
        self.ui.checkBox_VTA_From_E.stateChanged.connect(
            lambda: self.show_menu_item_on_checkbox_click(self.ui.checkBox_VTA_From_E,
                                                          self.ui.widget_Activation_Threshold_3))

        self.ui.lineEdit_Activation_Threshold.editingFinished.connect(
            lambda: self.check_lineedit_if_list_entered_absolute_float(self.ui.lineEdit_Activation_Threshold.text(),
                                                                       self.ui.lineEdit_Activation_Threshold))

        # unknown funcitons
        # self.ui.widget_Full_Field_IFFT.hide()
        self.ui.widget_Activation_Threshold_1.hide()
        self.ui.widget_Activation_Threshold_2.hide()
        self.ui.widget_Activation_Threshold_3.hide()

        # QPopUp
        self.filename = "{}/pop_up_control/dictionaries/dict_full_field_ifft.py".format(self.rel_folder)

        self.ui.pushButton_Save.clicked.connect(
            lambda: self.saveCloseWindow(self.output_dict(), self.filename))
        self.ui.pushButton_Cancel.clicked.connect(lambda: self.closeWindow())

    def output_dict(self):
        output_dict = {
            't_step_end': self.set_default_values(12000, self.ui.spinBox_T_Step_End.value()),
            'VTA_from_divE': self.corrector(0, self.ui.checkBox_VTA_From_DivE.checkState()),
            'VTA_from_NEURON': self.corrector(0, self.ui.checkBox_VTA_From_NEURON.checkState()),
            'VTA_from_E': self.corrector(0, self.ui.checkBox_VTA_From_E.checkState()),
            'Activation_threshold_VTA': self.set_default_values(0,
                                                                self.ui.lineEdit_Activation_Threshold.text())
        }
        return output_dict

    def saveDict(self):
        self.saveCloseWindow(self.output_dict(), self.filename)
