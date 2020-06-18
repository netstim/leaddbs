import sys
from PyQt5.QtWidgets import QApplication, QWidget
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5 import *

blue = "#309AFE"
green = "#89C843"
red = "#FD0100"


class qt_gui(QWidget):
    def __init__(self):
        super().__init__()
        # self.setFixedSize(1600, 600)
        # Create image and paint white
        # self.image = QImage(1600, 600, QImage.Format_RGB32)
        # self.image.fill(Qt.white)
        # self.update()

        self.canvas = QtGui.QPixmap(1600, 600)
        self.canvas.fill(Qt.white)
        self.drawSomething()

    def drawSomething(self):

        ref_y_pos = 600/ 2
        ref_x_pos = 50
        size_x = 200
        size_y = 100
        shift_x = 100
        shift_y = 50

        # define nodes
        mc_n = [[200, 225], [250, 270], [250, 300], [250, 330], [200, 375]]
        th_n = [[550, 125], [550, 170], [350, 170], [350, 125]]
        str_n = [[550, 420], [550, 450], [550, 480], [350, 480], [350, 420]]
        gpe_n = [[850, 270], [850, 300], [850, 330], [800, 350], [750, 350], [700, 350]]
        stn_n = [[1150, 280], [1150, 320], [1100, 350], [1050, 350], [1000, 350], [950, 330], [950, 300]]
        gpi_n = [[1250, 340], [1250, 320], [1250, 280], [1250, 260], [1295, 250]]
        sn_n = [[1450, 450], [1400, 500], [1300, 500], [1250, 480], [1250, 450]]

        # Set defaults
        if True:
            # h_mc-stn-gpi
            self.drawLine(stn_n[0], gpi_n[2], 'h', '', 0)
            self.drawTri(gpi_n[2], 'h', 'r', 0)
            # h_mc-stn-sn
            self.drawLine(stn_n[2], sn_n[2], 'h', '', 0, [stn_n[2][0], sn_n[2][1] + 40],
                          [sn_n[2][0], sn_n[2][1] + 40],
                          [sn_n[2][0], sn_n[2][1]])
            self.drawTri(sn_n[2], 'h', 'u', 0)
            # h_mc-stn
            self.drawLine(mc_n[2], stn_n[6], 'h', '', 0, [600, mc_n[2][1]], [600, mc_n[2][1] + 100],
                          [stn_n[4][0], mc_n[2][1] + 100], [stn_n[4][0], stn_n[4][1]])
            self.drawTri(stn_n[4], 'h', 'u', 0)

            # "i": ['MC_Str_GPe_GPi', 'MC_Str_GPe_STN_GPi', 'MC_Str_GPe_STN_SN']
            # i_gpe-stn
            self.drawLine(gpe_n[2], stn_n[5], 'i', '', 0)
            self.drawCircle(stn_n[5], "i", 0)
            # i_gpe-sn
            self.drawLine(gpe_n[3], sn_n[4], 'i', 'vh', 0)
            self.drawCircle(sn_n[4], "i", 0)
            # i_stn-gpi
            self.drawLine(stn_n[1], gpi_n[1], 'i', '', 0)
            self.drawTri(gpi_n[1], 'i', 'r', 0)
            # i_stn-sn
            # i_str-gpe
            self.drawLine(str_n[1], gpe_n[4], 'i', 'hv', 0)
            self.drawLine(stn_n[3], sn_n[1], 'i', '', 0, [stn_n[3][0], sn_n[1][1] + 60],
                          [sn_n[1][0], sn_n[1][1] + 60],
                          [sn_n[1][0], sn_n[1][1]])
            self.drawTri(sn_n[1], 'i', 'u', 0)
            self.drawLine(mc_n[4], str_n[3], 'i', 'vh', 0)
            self.drawTri(str_n[3], 'i', 'r', 0)
            self.drawLine(gpe_n[0], gpi_n[3], 'i', '', 0, [gpe_n[0][0] + 50, gpe_n[0][1]], [gpe_n[0][0] + 50, gpe_n[0][1] - 50], [gpi_n[3][0] - 50, gpe_n[0][1] - 50], [gpi_n[3][0] - 50, gpi_n[3][1]], [gpi_n[3][0], gpi_n[3][1]])
            self.drawCircle(gpi_n[3], "i", 0)

            # draw circles
            self.drawCircle(gpe_n[4], "i", 0)

            # "d": ['MC_Str_GPi', 'MC_Str_SN']
            # d_str-sn
            self.drawLine(str_n[2], sn_n[3], 'd', '', 0)
            self.drawCircle(sn_n[3], "d", 0)
            # d_mc-str
            self.drawLine(mc_n[3], str_n[4], 'd', '', 0, [mc_n[3][0] + 50, mc_n[3][1]],
                          [mc_n[3][0] + 50, str_n[4][1]],
                          [str_n[4][0], str_n[4][1]])
            self.drawTri(str_n[4], 'd', 'r', 0)
            self.drawLine(str_n[2], gpi_n[0], 'd', '', 0, [gpi_n[0][0] - 50, str_n[2][1]],
                          [gpi_n[0][0] - 50, gpi_n[0][1]], [gpi_n[0][0], gpi_n[0][1]])

            # draw circle
            self.drawCircle(gpi_n[0], "d", 0)

            # "e": ['STN_GPe', 'Th_MC'
            # e_mc-th
            self.drawLine(mc_n[0], th_n[3], 'e', 'vh', 0)
            self.drawTri(mc_n[0], 'e', 'b', 0)
            # e_th-sn
            # e_stn-gpe
            self.drawLine(gpe_n[1], stn_n[6], 'e', '', 0)
            self.drawTri(gpe_n[1], 'e', 'l', 0)
            self.drawLine(mc_n[1], th_n[2], 'e', '', 0, [mc_n[1][0] + 50, mc_n[1][1]],
                          [mc_n[1][0] + 50, th_n[2][1]],
                          [th_n[2][0], th_n[2][1]])
            self.drawTri(th_n[2], 'e', 'r', 0)
            # ib_mc-th
            # ib_th-gpi
            self.drawLine(th_n[0], gpi_n[4], 'ib', 'hv', 0)
            self.drawCircle(th_n[0], "ib", 0)
            self.drawCircle(th_n[1], "ib", 0)
            # ib_str-gpe
            self.drawLine(str_n[0], gpe_n[5], 'ib', 'hv', 0)
            self.drawCircle(str_n[0], "ib", 0)
            self.drawLine(th_n[1], sn_n[0], 'e', '', 0, [sn_n[0][0] + 50, th_n[1][1]],
                          [sn_n[0][0] + 50, sn_n[0][1]],
                          [sn_n[0][0], sn_n[0][1]])

        for dict_key in dict_connections:
            # draw lines
            # groups
            # "h": ['MC_STN', 'MC_STN_GPi', 'MC_STN_SN']
            result = [-1, -1]
            result[0] = dict_key
            result[1] = dict_connections[dict_key]

            if 'hdp' in result[0].lower():
                if 'stn_gpi' in result[0].lower():
                    # h_mc-stn-gpi
                    self.drawLine(stn_n[0], gpi_n[2], 'h', '', result[1])
                    self.drawTri(gpi_n[2], 'h', 'r', result[1])
                if 'stn_sn' in result[0].lower():
                    # h_mc-stn-sn
                    self.drawLine(stn_n[2], sn_n[2], 'h', '', result[1], [stn_n[2][0], sn_n[2][1] + 40], [sn_n[2][0], sn_n[2][1] + 40],
                          [sn_n[2][0], sn_n[2][1]])
                    self.drawTri(sn_n[2], 'h', 'u', result[1])
                if 'mc_stn' in result[0].lower():
                    # h_mc-stn
                    self.drawLine(mc_n[2], stn_n[6], 'h', '', result[1], [600, mc_n[2][1]], [600, mc_n[2][1] + 100],
                              [stn_n[4][0], mc_n[2][1] + 100], [stn_n[4][0], stn_n[4][1]])
                    self.drawTri(stn_n[4], 'h', 'u', result[1])

            # "i": ['MC_Str_GPe_GPi', 'MC_Str_GPe_STN_GPi', 'MC_Str_GPe_STN_SN']
            elif 'indirect' in result[0].lower():
                if 'gpe_stn' in result[0].lower():
                    # i_gpe-stn
                    self.drawLine(gpe_n[2], stn_n[5], 'i', '', result[1])
                    self.drawCircle(stn_n[5], "i", result[1])
                if 'gpe_sn' in result[0].lower():
                    # i_gpe-sn
                    self.drawLine(gpe_n[3], sn_n[4], 'i', 'vh', result[1])
                    self.drawCircle(sn_n[4], "i", result[1])
                if 'stn_gpi' in result[0].lower():
                    # i_stn-gpi
                    self.drawLine(stn_n[1], gpi_n[1], 'i', '', result[1])
                    self.drawTri(gpi_n[1], 'i', 'r', result[1])
                # i_stn-sn
                if 'str_gpe' in result[0].lower():
                    # i_str-gpe
                    self.drawLine(str_n[1], gpe_n[4], 'i', 'hv', result[1])
                if 'stn_sn' in result[0].lower():
                    self.drawLine(stn_n[3], sn_n[1], 'i', '', result[1], [stn_n[3][0], sn_n[1][1] + 60], [sn_n[1][0], sn_n[1][1] + 60],
                              [sn_n[1][0], sn_n[1][1]])
                    self.drawTri(sn_n[1], 'i', 'u', result[1])
                if 'mc_str' in result[0].lower():
                    self.drawLine(mc_n[4], str_n[3], 'i', 'vh', result[1])
                    self.drawTri(str_n[3], 'i', 'r', result[1])
                if 'gpe_gpi' in result[0].lower():
                    self.drawLine(gpe_n[0], gpi_n[3], 'i', '', result[1], [gpe_n[0][0] + 50, gpe_n[0][1]],
                              [gpe_n[0][0] + 50, gpe_n[0][1] - 50], [gpi_n[3][0] - 50, gpe_n[0][1] - 50],
                              [gpi_n[3][0] - 50, gpi_n[3][1]], [gpi_n[3][0], gpi_n[3][1]])
                    self.drawCircle(gpi_n[3], "i", result[1])

                    # draw circles
                    self.drawCircle(gpe_n[4], "i", result[1])

            # "d": ['MC_Str_GPi', 'MC_Str_SN']
            elif 'direct' in result[0].lower():
                if 'str_sn' in result[0].lower():
                    # d_str-sn
                    self.drawLine(str_n[2], sn_n[3], 'd', '', result[1])
                    self.drawCircle(sn_n[3], "d", result[1])
                if 'mc_str' in result[0].lower():
                    # d_mc-str
                    self.drawLine(mc_n[3], str_n[4], 'd', '', result[1], [mc_n[3][0] + 50, mc_n[3][1]], [mc_n[3][0] + 50, str_n[4][1]],
                              [str_n[4][0], str_n[4][1]])
                    self.drawTri(str_n[4], 'd', 'r', result[1])
                if 'str_gpi' in result[0].lower():
                    self.drawLine(str_n[2], gpi_n[0], 'd', '', result[1], [gpi_n[0][0] - 50, str_n[2][1]], [gpi_n[0][0] - 50, gpi_n[0][1]], [gpi_n[0][0], gpi_n[0][1]])

                    # draw circle
                    self.drawCircle(gpi_n[0], "d", result[1])

                # "e": ['STN_GPe', 'Th_MC'
            elif 'excitatory' in result[0].lower():
                if 'th_mc' in result[0].lower():
                    # e_mc-th
                    self.drawLine(mc_n[0], th_n[3], 'e', 'vh', result[1])
                    self.drawTri(mc_n[0], 'e', 'b', result[1])
                # e_th-sn
                if 'stn_gpe' in result[0].lower():
                    # e_stn-gpe
                    self.drawLine(gpe_n[1], stn_n[6], 'e', '', result[1])
                    self.drawTri(gpe_n[1], 'e', 'l', result[1])
                if 'mc_th' in result[0].lower():
                    self.drawLine(mc_n[1], th_n[2], 'e', '', result[1], [mc_n[1][0] + 50, mc_n[1][1]], [mc_n[1][0] + 50, th_n[2][1]],
                                  [th_n[2][0], th_n[2][1]])
                    self.drawTri(th_n[2], 'e', 'r', result[1])
            elif 'inhibitory' in result[0].lower():
                # ib_mc-th
                if 'gpi_th' in result[0].lower():
                    # ib_th-gpi
                    self.drawLine(th_n[0], gpi_n[4], 'ib', 'hv', result[1])
                    self.drawCircle(th_n[0], "ib", result[1])
                    self.drawCircle(th_n[1], "ib", result[1])
                if 'gpe_str' in result[0].lower():
                    # ib_str-gpe
                    self.drawLine(str_n[0], gpe_n[5], 'ib', 'hv', result[1])
                    self.drawCircle(str_n[0], "ib", result[1])
                if 'sn_th' in result[0].lower():
                    self.drawLine(th_n[1], sn_n[0], 'e', '', result[1], [sn_n[0][0] + 50, th_n[1][1]], [sn_n[0][0] + 50, sn_n[0][1]],
                                  [sn_n[0][0], sn_n[0][1]])

        # draw rectangles and text
        self.drawRectText([ref_x_pos, ref_y_pos], [size_x, size_y + 50], " Motor\nCortex")
        self.update()
        self.drawRectText([ref_x_pos + size_x + shift_x, ref_y_pos - (size_y + shift_y)], [size_x, size_y], "Thalamus")
        self.drawRectText([ref_x_pos + size_x + shift_x, ref_y_pos + (size_y + shift_y)], [size_x, size_y], "Striatum")
        self.drawRectText([ref_x_pos + 2 * (size_x + shift_x), ref_y_pos + 0 * (size_y + shift_y)], [size_x, size_y],
                          "GPe")
        self.drawRectText([ref_x_pos + 3 * (size_x + shift_x), ref_y_pos + 0 * (size_y + shift_y)], [size_x, size_y],
                          "STN")
        self.drawRectText([ref_x_pos + 4 * (size_x + shift_x), ref_y_pos + 0 * (size_y + shift_y)], [size_x, size_y],
                          "GPi")
        self.drawRectText([ref_x_pos + 4 * (size_x + shift_x), ref_y_pos + (size_y + shift_y)], [size_x, size_y], "SN")


        for dict_key in dict_connections:
            result=[-1,-1]
            result[0]=dict_key
            result[1]=dict_connections[dict_key]*4

            # draw lines
            # groups[]
            # "h": ['MC_STN', 'MC_STN_GPi', 'MC_STN_SN']
            #print(result)
            #print(result[0])
            if 'hdp' in result[0].lower():
                if 'stn_gpi' in result[0].lower():
                    # h_mc-stn-gpi
                    # self.drawLine(stn_n[0], gpi_n[2], 'h', '', result[1])
                    self.drawTri(gpi_n[2], 'h', 'r', result[1])
                if 'stn_sn' in result[0].lower():
                    # h_mc-stn-sn
                    # self.drawLine(stn_n[2], sn_n[2], 'h', '', result[1], [stn_n[2][0], sn_n[2][1] + 40], [sn_n[2][0], sn_n[2][1] + 40],
                    #       [sn_n[2][0], sn_n[2][1]])
                    self.drawTri(sn_n[2], 'h', 'u', result[1])
                if 'mc_stn' in result[0].lower():
                    # h_mc-stn
                    # self.drawLine(mc_n[2], stn_n[6], 'h', '', result[1], [600, mc_n[2][1]], [600, mc_n[2][1] + 100],
                    #           [stn_n[4][0], mc_n[2][1] + 100], [stn_n[4][0], stn_n[4][1]])
                    self.drawTri(stn_n[4], 'h', 'u', result[1])

            # "i": ['MC_Str_GPe_GPi', 'MC_Str_GPe_STN_GPi', 'MC_Str_GPe_STN_SN']
            elif 'indirect' in result[0].lower():
                if 'gpe_stn' in result[0].lower():
                    # i_gpe-stn
                    # self.drawLine(gpe_n[2], stn_n[5], 'i', '', result[1])
                    self.drawCircle(stn_n[5], "i", result[1])
                if 'gpe_sn' in result[0].lower():
                    # i_gpe-sn
                    # self.drawLine(gpe_n[3], sn_n[4], 'i', 'vh', result[1])
                    self.drawCircle(sn_n[4], "i", result[1])
                if 'stn_gpi' in result[0].lower():
                    # i_stn-gpi
                    # self.drawLine(stn_n[1], gpi_n[1], 'i', '', result[1])
                    self.drawTri(gpi_n[1], 'i', 'r', result[1])
                # i_stn-sn
                # if 'str_gpe' in result[0].lower():
                    # i_str-gpe
                    # self.drawLine(str_n[1], gpe_n[4], 'i', 'hv', result[1])
                if 'stn_sn' in result[0].lower():
                    # self.drawLine(stn_n[3], sn_n[1], 'i', '', result[1], [stn_n[3][0], sn_n[1][1] + 60], [sn_n[1][0], sn_n[1][1] + 60], [sn_n[1][0], sn_n[1][1]])
                    self.drawTri(sn_n[1], 'i', 'u', result[1])
                if 'mc_str' in result[0].lower():
                    # self.drawLine(mc_n[4], str_n[3], 'i', 'vh', result[1])
                    self.drawTri(str_n[3], 'i', 'r', result[1])
                if 'gpe_gpi' in result[0].lower():
                    # self.drawLine(gpe_n[0], gpi_n[3], 'i', '', result[1], [gpe_n[0][0] + 50, gpe_n[0][1]],
                    #           [gpe_n[0][0] + 50, gpe_n[0][1] - 50], [gpi_n[3][0] - 50, gpe_n[0][1] - 50],
                    #           [gpi_n[3][0] - 50, gpi_n[3][1]], [gpi_n[3][0], gpi_n[3][1]])
                    self.drawCircle(gpi_n[3], "i", result[1])
                    # draw circles
                    self.drawCircle(gpe_n[4], "i", result[1])

            # "d": ['MC_Str_GPi', 'MC_Str_SN']
            elif 'direct' in result[0].lower():
                if 'str_sn' in result[0].lower():
                    # d_str-sn
                    # self.drawLine(str_n[2], sn_n[3], 'd', '', result[1])
                    self.drawCircle(sn_n[3], "d", result[1])
                if 'mc_str' in result[0].lower():
                    # d_mc-str
                    # self.drawLine(mc_n[3], str_n[4], 'd', '', result[1], [mc_n[3][0] + 50, mc_n[3][1]], [mc_n[3][0] + 50, str_n[4][1]],
                    #           [str_n[4][0], str_n[4][1]])
                    self.drawTri(str_n[4], 'd', 'r', result[1])
                if 'str_gpi' in result[0].lower():
                    # self.drawLine(str_n[2], gpi_n[0], 'd', '', result[1], [gpi_n[0][0] - 50, str_n[2][1]], [gpi_n[0][0] - 50, gpi_n[0][1]], [gpi_n[0][0], gpi_n[0][1]])

                    # draw circle
                    self.drawCircle(gpi_n[0], "d", result[1])

                # "e": ['STN_GPe', 'Th_MC'
            elif 'excitatory' in result[0].lower():
                if 'th_mc' in result[0].lower():
                    # e_mc-th
                    # self.drawLine(mc_n[0], th_n[3], 'e', 'vh', result[1])
                    self.drawTri(mc_n[0], 'e', 'b', result[1])
                # e_th-sn
                if 'stn_gpe' in result[0].lower():
                    # e_stn-gpe
                    # self.drawLine(gpe_n[1], stn_n[6], 'e', '', result[1])
                    self.drawTri(gpe_n[1], 'e', 'l', result[1])
                if 'mc_th' in result[0].lower():
                    # self.drawLine(mc_n[1], th_n[2], 'e', '', result[1], [mc_n[1][0] + 50, mc_n[1][1]], [mc_n[1][0] + 50, th_n[2][1]],
                    #               [th_n[2][0], th_n[2][1]])
                    self.drawTri(th_n[2], 'e', 'r', result[1])
            elif 'inhibitory' in result[0].lower():
                # ib_mc-th
                if 'gpi_th' in result[0].lower():
                    # ib_th-gpi
                    # self.drawLine(th_n[0], gpi_n[4], 'ib', 'hv', result[1])
                    self.drawCircle(th_n[0], "ib", result[1])
                    self.drawCircle(th_n[1], "ib", result[1])
                if 'gpe_str' in result[0].lower():
                    # ib_str-gpe
                    # self.drawLine(str_n[0], gpe_n[5], 'ib', 'hv', result[1])
                    self.drawCircle(str_n[0], "ib", result[1])
                # if 'sn_th' in result[0].lower():
                    # self.drawLine(th_n[1], sn_n[0], 'e', '', result[1], [sn_n[0][0] + 50, th_n[1][1]], [sn_n[0][0] + 50, sn_n[0][1]],
                    #               [sn_n[0][0], sn_n[0][1]])

        # housekeeping
        size = 0.7
        self.housekeep([600, 50], "excitatory")
        self.drawLine([610, 100], [700, 100], 'e', '', size)
        self.drawTri([710, 100], 'e', 'r', size)
        self.housekeep([800, 50], "inhibitory")
        self.drawLine([810, 100], [900, 100], 'ib', '', size)
        self.drawCircle([900, 100], 'ib', size)
        self.housekeep([1020, 50], "HDP")
        self.drawLine([1010, 100], [1100, 100], 'h', '', size)
        self.housekeep([1200, 50], "indirect")
        self.drawLine([1200, 100], [1300, 100], 'i', '', size)
        self.housekeep([1410, 50], "direct")
        self.drawLine([1400, 100], [1500, 100], 'd', '', size)

        # Draw activation
        if 'stn' in electrodeLoc.lower():
            self.drawBolt([stn_n[0][0]-10, stn_n[0][1]-30])
        elif 'th' in electrodeLoc.lower():
            self.drawBolt([th_n[0][0]-10, th_n[0][1]-20])
        elif 'gpi' in electrodeLoc.lower():
            self.drawBolt([gpi_n[0][0]+190, gpi_n[0][1]-90])
        elif 'sn' in electrodeLoc.lower():
            self.drawBolt([sn_n[0][0]-10, sn_n[0][1]-50])
        elif 'gpe' in electrodeLoc.lower():
            self.drawBolt([gpe_n[0][0]-10, gpe_n[0][1]-20])

    def drawRectText(self, p, s, t):

        # hire a painter
        painter = QPainter(self.canvas)

        # buy a pen
        pen = QPen()
        pen.setWidth(5)

        # give pen to painter
        painter.setPen(pen)

        # draw rectangle
        painter.fillRect(QRect(p[0], p[1] - 0.5 * s[1], s[0], s[1]), QColor("White"))
        painter.drawRoundedRect(QRect(p[0], p[1] - 0.5 * s[1], s[0], s[1]), 10, 10)

        # write text
        font_size = 30
        painter.setFont(QFont("Gothic", font_size))

        text_len = len(t)
        if len(t) > 8:
            text_len = 5

        painter.drawText(
            QRect(p[0] + s[0] / 2 - (text_len * font_size) / 2.7, p[1] - 0.5 * (s[1] - font_size), s[0], s[1]), 0, t)

        painter.end()

    def drawLine(self, start, end, pathway, orient, size, *points):
        # house keeping
        color = "Black"
        if pathway == 'd':
            color = red
        elif pathway == 'i':
            color = green
        elif pathway == 'ib' or pathway == 'e':
            color = "Black"
        elif pathway == 'h':
            color = blue

        # hire a painter
        painter = QPainter(self.canvas)

        # buy a pen
        pen = QPen(Qt.green, size*10, Qt.SolidLine, Qt.RoundCap, Qt.RoundJoin)
        pen.setColor(QColor(color))

        # give pen to painter
        painter.setPen(pen)

        # draw path
        drawingPath = QPainterPath()
        if points:
            drawingPath.moveTo(start[0], start[1])
            for point in points:
                drawingPath.lineTo(point[0], point[1])
        else:
            if orient == 'vh':
                drawingPath.moveTo(start[0], start[1])
                drawingPath.lineTo(start[0], end[1])
                drawingPath.lineTo(end[0], end[1])
            elif orient == 'hv':
                drawingPath.moveTo(start[0], start[1])
                drawingPath.lineTo(end[0], start[1])
                drawingPath.lineTo(end[0], end[1])
            else:
                drawingPath.moveTo(start[0], start[1])
                drawingPath.lineTo(start[0], end[1])
                drawingPath.lineTo(end[0], end[1])

        painter.drawPath(drawingPath)

        painter.end()

    def drawCircle(self, pos, pathway, size):
        # house keeping
        color = "Black"
        if pathway == 'd':
            color = red
        elif pathway == 'i':
            color = green
        elif pathway == 'ib' or pathway == 'e':
            color = "Black"
        elif pathway == 'h':
            color = blue

        # hire a painter
        painter = QPainter(self.canvas)

        # buy a pen
        pen = QPen(Qt.green, size*10, Qt.SolidLine, Qt.SquareCap, Qt.RoundJoin)
        pen.setColor(QColor(color))

        # give pen to painter
        painter.setPen(pen)

        painter.setBrush(QBrush(QColor(color), Qt.SolidPattern))
        painter.drawEllipse(pos[0] - 5, pos[1] - 5, 10, 10)

        painter.end()
        self.update()

    def housekeep(self, pos, t):
        # hire a painter
        painter = QPainter(self.canvas)

        # buy a pen
        pen = QPen(Qt.black, 5, Qt.SolidLine, Qt.SquareCap, Qt.RoundJoin)
        # pen.setColor(QColor(color))

        # give pen to painter
        painter.setPen(pen)

        # write text
        font_size = 20
        painter.setFont(QFont("Gothic", font_size))
        painter.drawText(QRect(pos[0], pos[1], 30 * len(t), 40), 0, t)


        painter.end()
        self.update()

    def drawTri(self, pos, pathway, orient, size):
        # house keeping
        color = "Black"
        if pathway == 'd':
            color = red
        elif pathway == 'i':
            color = green
        elif pathway == 'ib' or pathway == 'e':
            color = "Black"
        elif pathway == 'h':
            color = blue

        # hire a painter
        painter = QPainter(self.canvas)

        # buy a pen
        pen = QPen(Qt.black, size*10, Qt.SolidLine, Qt.SquareCap, Qt.RoundJoin)
        # pen.setColor(QColor(color))

        # give pen to painter
        painter.setPen(pen)

        # draw triangle
        i = (size+1)*10
        # draw path
        drawingPath = QPainterPath()
        if orient == "u":
            drawingPath.moveTo(pos[0], pos[1])
            drawingPath.lineTo(pos[0] + i / 2, pos[1] + i)
            drawingPath.lineTo(pos[0] - i / 2, pos[1] + i)
        elif orient == "b":
            drawingPath.moveTo(pos[0], pos[1])
            drawingPath.lineTo(pos[0] + i / 2, pos[1] - i)
            drawingPath.lineTo(pos[0] - i / 2, pos[1] - i)
        elif orient == "l":
            drawingPath.moveTo(pos[0], pos[1])
            drawingPath.lineTo(pos[0] + i, pos[1] - i / 2)
            drawingPath.lineTo(pos[0] + i, pos[1] + i / 2)
        elif orient == "r":
            drawingPath.moveTo(pos[0], pos[1])
            drawingPath.lineTo(pos[0] - i, pos[1] - i / 2)
            drawingPath.lineTo(pos[0] - i, pos[1] + i / 2)

        painter.fillPath(drawingPath, QBrush(QColor(color)))

        painter.end()
        self.update()

    def drawBolt(self, pos):
        # hire a painter
        painter = QPainter(self.canvas)

        # buy a pen
        pen = QPen(Qt.black, 10, Qt.SolidLine, Qt.SquareCap, Qt.RoundJoin)
        # pen.setColor(QColor(color))

        # give pen to painter
        painter.setPen(pen)

        # draw bolt
        drawingPath = QPainterPath()
        drawingPath.moveTo(pos[0], pos[1])
        drawingPath.lineTo(pos[0]-50, pos[1]+60)
        drawingPath.lineTo(pos[0]-30, pos[1]+60)
        drawingPath.lineTo(pos[0]-50, pos[1]+110)
        drawingPath.lineTo(pos[0], pos[1]+50)
        drawingPath.lineTo(pos[0]-20, pos[1]+50)


        drawingPath.lineTo(pos[0], pos[1])

        painter.fillPath(drawingPath, QBrush(Qt.yellow))

        painter.end()
        self.update()

    def saveImage(self, fileName, fileFormat):
        self.canvas.save(fileName, fileFormat)


class Button(QWidget):
    def __init__(self):
        super().__init__()
        # Set layout
        layout = QtWidgets.QVBoxLayout()

        # Create label to hold QPixmap
        canvas = qt_gui()
        label = QtWidgets.QLabel()
        label.setPixmap(canvas.canvas)

        # Create instance of paint canvas and add to layout
        layout.addWidget(label)

        # Create new horizontal layout for buttons and add buttons to layout
        hlayout = QtWidgets.QHBoxLayout()

        # Create spacer
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        hlayout.addItem(spacerItem)

        self.savebutton = QtWidgets.QPushButton()
        self.savebutton.setText("Save")
        hlayout.addWidget(self.savebutton)

        self.closebutton = QtWidgets.QPushButton()
        self.closebutton.setText("Close")
        hlayout.addWidget(self.closebutton)

        # Add horizontal layout to vertical layout
        layout.addLayout(hlayout)

        self.setLayout(layout)

        # Call to events function
        self.events()

    def events(self):
        self.savebutton.clicked.connect(lambda: self.save_image())
        self.closebutton.clicked.connect(lambda: sys.exit())

    def save_image(self):
        print("this is working")
        paintCanvas = qt_gui()
        paintCanvas.saveImage("trial.png", "PNG")


def main():
    app = QApplication(sys.argv)
    button = Button()
    button.show()
    sys.exit(app.exec_())


def drawResults(resultsInput, activation):
    #print(resultsInput)
    global dict_connections, electrodeLoc
    electrodeLoc = activation
    dict_connections = resultsInput
    main()
