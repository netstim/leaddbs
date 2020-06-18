from PyQt5.QtWidgets import QWidget
from functionalities import Functionalities
from pop_up_ui import pop_up_CPE_Active


class PopUpCPEActive(Functionalities):
    def __init__(self, parent):
        self.main_win = QWidget()
        self.ui = pop_up_CPE_Active.Ui_centralwidget()
        self.ui.setupUi(self.main_win)

        # Create main window object
        self.mainWindow = parent


        # mini dictionary
        self.filename = "{}/pop_up_control/dictionaries/dict_cpe_active.py".format(self.rel_folder())

        # Save and cancel
        self.ui.pushButton_Save.clicked.connect(
            lambda: self.saveCloseWindow(self.output_dict(), self.filename))
        self.ui.pushButton_Cancel.clicked.connect(lambda: self.closeWindow())

    def output_dict(self):
        output_dict = {
            'beta': self.set_default_values(0, self.ui.doubleSpinBox_Alpha.value()),
            'K_A': self.set_default_values(0, self.ui.doubleSpinBox_K_S.value()),
            'beta_ground': self.set_default_values(0, self.ui.doubleSpinBox_Alpha_Ground.value()),
            'K_A_ground': self.set_default_values(0, self.ui.doubleSpinBox_K_S_Ground.value()),
        }

        return output_dict

    def saveDict(self):
        self.saveCloseWindow(self.output_dict(), self.filename)
