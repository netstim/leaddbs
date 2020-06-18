import ast
import os
import re
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

import Header


class Functionalities:
    def rel_folder(self):
        self.rel_folder = 'GUI_tree_files' # the relative path could be changed from here at anytime
        return self.rel_folder

    # QPopUp functionality
    def pop_up_push_button(self, push_button, window, checkbox):
        if push_button.text() == "Save":
            window.close()
        elif push_button.text() == "Cancel":
            window.close()
            if checkbox.checkState() == 0:
                checkbox.setCheckState(2)
            else:
                checkbox.setCheckState(0)

    def pop_up(self, check_box_state, pop_up_window, main_window):
        if check_box_state.checkState() == 2:
            pop_up_window.showWindow()
            # main_window.setDisabled(True)
        else:
            pass

    def pop_up_reverse(self, check_box_state, pop_up_window, main_window):
        if check_box_state.checkState() == 0:
            pop_up_window.showWindow()
            # main_window.setDisabled(True)
        else:
            pass

    def return_main_window_function(self, button, pop_up_window, main_window):
        if button.text == "Save":
            pop_up_window.closeWindow()
            # main_window.setDisabled(False)
        elif button.text == "Cancel":
            pop_up_window.closeWindow()
            # main_window.setDisabled(False)

    def show_menu_item_on_checkbox_click(self, widget, menu):
        if widget.checkState() == 2:
            menu.show()
        else:
            menu.hide()

    def hide_menu_item_on_checkbox_click(self, widget, menu):
        if widget.checkState() == 2:
            menu.hide()
        else:
            menu.show()

    def show_menu_item_on_checkbok_click_double_menu(self, widget, menu1, menu2):
        if widget.checkState() == 2:
            menu1.show()
            menu2.hide()
        else:
            menu1.hide()
            menu2.show()

    def enable_disable_menu_item_on_checkbok_click_double_menu(self, widget, menu1, menu2):
        if widget.checkState() == 2:
            menu1.show()
            menu2.hide()
        elif widget.checkState() == 0:
            menu1.hide()
            menu2.show()

    def spectrum_truncation_method_output(self, text, widget, out1, out2):
        if widget.currentText() == text:
            return out1.text()
        else:
            return out2.value()

    def show_menu_item_on_combobox_state_change(self, text, widget, menu):
        if widget.currentText() == text:
            menu.hide()
        else:
            menu.show()

    def hide_menu_item_on_combobox_state_change(self, text, widget, menu):
        if widget.currentText() == text:
            menu.show()
        else:
            menu.hide()

    def show_menu_item_lineedit(self, x, item):
        if x.text() == "":
            print(x.text())
            item.hide()
        else:
            item.show()

    # Help Message Box Functions
    def info(self, title, info):
        msg_box = QMessageBox()
        msg_box.setWindowTitle(title)
        msg_box.setText(info)
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setStandardButtons(QMessageBox.Ok)
        # msg_box.setDefaultButton(QMessageBox.Ok)
        # msg_box.setInformativeText("Informative text")
        # msg_box.setDetailedText("This text is initially hidden")
        x = msg_box.exec_()

    # CheckBoxes/Spinners writing to dictionary functions
    def corrector(self, default, arg):
        if arg != "":
            if arg > 1:
                arg -= 1
                return arg
            else:
                return arg
        else:
            return default

    def anti_corrector(self, arg):
        if arg == 1:
            arg += 1
            return arg
        else:
            return arg

    def convert_to_boolean(self, default, arg):
        if arg != "":
            if arg >= 1:
                return True
            else:
                return False
        else:
            return default

    def set_default_values(self, default, value):
        if value == "":
            return default
        else:
            return value

    def set_default_text(self, default, value):
        if value == "":
            return "'" + default + "'"
        else:
            return "'" + value + "'"

    def connected_check_boxes(self, checkbox1, checkbox2):
        if checkbox1.checkState() == 2:
            checkbox2.setCheckState(0)

    def connected_check_boxes_3(self, checkbox1, checkbox2, checkbox3):
        if checkbox1.checkState() == 2:
            checkbox2.setCheckState(0)
            checkbox3.setCheckState(0)

    def disable_widget(self, control, widget):
        if control.checkState() == 2:
            widget.setDisabled(True)
            widget.setCheckState(0)
        else:
            widget.setDisabled(False)

    def check_lineedit_if_list_entered(self, entry, lineedit_widget):
        alphabets = "abcdfghijklmpqrstuvwxyz"
        for char in entry:
            if char in alphabets:
                lineedit_widget.setText('')
                break
        else:
            case_insensitive_none = re.compile(re.escape("none"), re.IGNORECASE)
            entry = case_insensitive_none.sub("None", entry)
            for char in entry:
                if char in '[]1234567890noe':
                    try:
                        entry = entry.replace(']', '').replace('[', '').replace(' ', ',').replace(',,', ',').strip(',')
                        entry_list = "[{}]".format(entry)
                        entry_f = ast.literal_eval(entry_list)
                        lineedit_widget.setText("{}".format(entry_f))
                        break
                    except:
                        lineedit_widget.setText('')
                        break
            else:
                lineedit_widget.setText("{}".format(entry))

    def check_lineedit_if_value_entered(self, entry, lineedit_widget):
        alphabets = "abcdefghijklmnopqrstuvwxyz"
        for char in entry:
            if char in alphabets:
                lineedit_widget.setText('')
                break
        else:
            for char in entry:
                if char in "[]" or char in ",":
                    lineedit_widget.setText('')
                    break
            else:
                lineedit_widget.setText("{}".format(entry))

    def check_lineedit_if_list_entered_absolute(self, entry, lineedit_widget):
        alphabets = "abcdfghijklmpqrstuvwxyz"
        for char in entry:
            if char in alphabets:
                lineedit_widget.setText('')
                break
        else:
            case_insensitive_none = re.compile(re.escape("none"), re.IGNORECASE)
            entry = case_insensitive_none.sub("None", entry)
            for char in entry:
                if char in '[]1234567890noe':
                    try:
                        entry = entry.replace(']', '').replace('[', '').replace(' ', ',').replace(',,', ',').replace('-', '').strip(',')
                        entry_list = "[{}]".format(entry)
                        entry_f = ast.literal_eval(entry_list)
                        lineedit_widget.setText("{}".format(entry_f))
                        break
                    except:
                        lineedit_widget.setText('')
                        break
            else:
                lineedit_widget.setText("{}".format(entry))

    def check_lineedit_if_list_entered_absolute_float(self, entry, lineedit_widget):
        alphabets = "abcdfghijklmpqrstuvwxyz"
        for char in entry:
            if char in alphabets:
                lineedit_widget.setText('')
                break
        else:
            case_insensitive_none = re.compile(re.escape("none"), re.IGNORECASE)
            entry = case_insensitive_none.sub("None", entry)
            for char in entry:
                if char in '[]1234567890noe':
                    try:
                        entry = entry.replace(']', '').replace('[', '').replace(' ', ',').replace(',,', ',').replace('-', '').strip(',')
                        entry_list = "[{}]".format(entry)
                        entry_f = ast.literal_eval(entry_list)
                        lineedit_widget.setText("{}".format(entry_f).strip('[]'))
                        break
                    except:
                        lineedit_widget.setText('')
                        break
            else:
                lineedit_widget.setText("{}".format(entry))

    def set_lineedit_to_list(self, entry):
        if isinstance(entry, float) or isinstance(entry, int):
            return "[{}]".format(entry)
        else:
            return entry

    def get_lineedit_entry(self, parameter, lineedit_widget):
        if isinstance(parameter, str):
            lineedit_widget.setText(parameter)
        else:
            parameter = "{}".format(parameter)
            lineedit_widget.setText(parameter)

    def get_entry_truncation_parameter(self, parameter, widget, object1, object2):
        if widget.currentText() == 'Octave Band Method':
            if isinstance(parameter, str):
                object1.setText(parameter)
            else:
                parameter = "{}".format(parameter)
                object1.setText(parameter)
        else:
            if isinstance(parameter, int):
                object2.setValue(parameter)
            else:
                parameter = int(parameter)
                object2.setValue(parameter)

    def set_current_file_name(self, current_file_name):
        self.current_file_name = current_file_name

    def py_check(self, filename):
        if filename.endswith(".py"):
            return filename
        else:
            return filename + ".py"

    def txt_check(self, filename):
        if filename.endswith(".txt"):
            return filename
        else:
            return filename + ".txt"

    def loadStylesheet(self, filename):
        file = QFile(filename)
        file.open(QFile.ReadOnly | QFile.Text)
        stylesheet = file.readAll()
        QApplication.instance().setStyleSheet(str(stylesheet, encoding="utf-8"))

    def showWindow(self):
        self.main_win.show()

    def saveCloseWindow(self, output_dict, filename):
        # print("Window was saved and closed")
        self.dict_write(output_dict, filename)
        try:
            self.main_win.close()
        except:
            pass

    def closeWindow(self):
        # print("Window was closed")
        self.main_win.close()

    def dict_write(self, output_dict, filename):
        # print("Functioanlities:: Start saving")
        with open(filename, 'w') as test_dict:
            test_dict.write(Header.header)
            test_dict.write("d = {\n")
            for x in output_dict:
                test_dict.write("    '{}': {},\n".format(x, output_dict[x]))
            test_dict.write("    }\n")

        # with open("GUI_tree_files/last_save.py", 'w') as last_save:
        #     last_save.write("d = {\n")
        #     for x in output_dict:
        #         last_save.write("    '{}': {},\n".format(x, output_dict[x]))
        #     last_save.write("    }\n")
        #
        # self.info("Run", "Setup state has been written to {}.".format(filename))

    def display(self, label, path1, path2, mode=1):
        if os.path.exists(path1):
            images_dir_name = os.path.join(os.getcwd(), path1)
            # label.setScaledContents(True)
            sizePolicy = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
            label.setSizePolicy(sizePolicy)
            label.setAlignment(Qt.AlignCenter)
            if mode == 2:
                label.setPixmap(QPixmap(images_dir_name))
            else:
                label.setPixmap(QPixmap(images_dir_name).scaled(1280, 640, Qt.KeepAspectRatio))
        elif os.path.exists(path2):
            images_dir_name = os.path.join(os.getcwd(), path2)
            # label.setScaledContents(True)
            sizePolicy = QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
            label.setSizePolicy(sizePolicy)
            label.setAlignment(Qt.AlignCenter)
            if mode == 2:
                label.setPixmap(QPixmap(images_dir_name))
            else:
                label.setPixmap(QPixmap(images_dir_name).scaled(1280, 640, Qt.KeepAspectRatio))
        else:
            pass
