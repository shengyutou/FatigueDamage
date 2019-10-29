#!/usr/bin/env python
# _*_coding:utf-8 _*_
"""
@Time    : 10/29/19 2:30 PM
@Author  : Louis Wang
@Email   : louis1990.wang@outlook.com
@FileName: SN_Cal_GUI.py
@Software: PyCharm
"""


import sys
import time
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
from math import log10, sqrt
from GUI import SNCurve

MainWindow_main = SNCurve.Ui_SN_Curve


class Main(QtWidgets.QMainWindow, MainWindow_main):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        MainWindow_main.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QtGui.QIcon('SN.png'))
        self.Cal_Plot.clicked.connect(self.sn_cal)

    def sn_cal(self):
        """
        Calculation of the SN curve
        :return:
        """
        # Para Read
        try:
            self.mat = str(self.l_mat.currentText())
            self.Rp = float(self.l_Rp.text())
            self.Rm = float(self.l_Rm.text())
            self.t = float(self.l_t.text())
            self.Rz = float(self.l_Rz.text())
            self.j = float(self.l_j.text())
            self.j_0 = float(self.l_j_0.text())
            self.S_pu = float(self.l_S_pu.text())
            self.gamma_M = float(self.l_gamma_M.text())
            self.sigm_b = self.Rm * 1.06
            self.R = -1
            if self.check_t_c.isChecked():
                self.sign_t = 1
            else:
                self.sign_t = 0
        except Exception as e:
            with open('ErrorLog.txt', 'a') as f:
                f.write(time.strftime("%Y-%m-%d-%H-%M", time.localtime(time.time())) + '\n')
                f.write(repr(e))
                f.close()
                QMessageBox.warning(self, 'Error', '参数输入错误，请检查！\n')
            return

        # Calculation
        if self.mat == 'GGG':
            self.M = 0.00035 * self.sigm_b + 0.08
        else:
            self.M = 0.00035 * self.sigm_b + 0.05

        # surface roughness factor
        F_o = 1 - 0.22 * (log10(self.Rz) ** 0.64) * log10(self.sigm_b) + 0.45 * (log10(self.Rz) ** 0.53)

        # notch factor
        alph_k = 1
        n = 1
        belt_k = alph_k / n

        # total influence factor fok
        F_ok = sqrt(belt_k ** 2 - 1 + 1 / (F_o ** 2))

        # fatigue strength of specimen
        if self.mat == 1:
            sigm_w = 0.27 * self.sigm_b + 100
        else:
            sigm_w = 0.27 * self.sigm_b + 85

        # fatigue strength of component
        sigm_wk = sigm_w / F_ok

        # slopes of SN curve m1 and m2
        self.m1 = 5.5 / (F_ok ** 2) + 6
        self.m2 = 2 * self.m1 - 1

        # factor for influence of mean stress
        u = 1 / (self.M + 1) * sigm_wk / self.sigm_b
        a = (1 + self.R) / (1 - self.R) * sigm_wk / self.sigm_b
        p = (1 / (self.M + 1) - 1 + u ** 2) / (u ** 2 - u)
        if a == 0:
            Fm = 1
        else:
            if p <= 1:
                Fm = -1 * (1 + p * a) / (2 * a ** 2 * (1 - p)) + sqrt(
                    1 / (1 - p) / a ** 2 + ((1 + p * a) / 2 / a ** 2 / (1 - p)) ** 2)
            else:
                Fm = -1 * (1 + p * a) / (2 * a ** 2 * (1 - p)) - sqrt(
                    1 / (1 - p) / a ** 2 + ((1 + p * a) / 2 / a ** 2 / (1 - p)) ** 2)

        # stress amplitude at knee of SN curve
        sigm_A = sigm_wk * Fm

        # number of load cycles at knee of SN curve
        self.N_D = 10 ** (6.8 - 3.6 * (1 / self.m1))

        # upgrading factors
        # quelity level
        S_d = 0.85 ** (self.j - self.j_0)
        # thickness-dependent tensile value Rm
        S_t = (self.t / 25) ** ((-0.15) * self.sign_t)
        # total upgrading factor
        S = self.S_pu * S_d * S_t

        # upgraded stress amplitude at knee of SN curve
        self.sigm_d = sigm_A * S / self.gamma_M
        # upper limit of fatigue life line
        self.sigm_1 = self.Rp * (1 - self.R) / self.gamma_M
        # number of load cycles at upper fatigue limit
        self.N_1 = self.N_D * (2 * self.sigm_d / self.sigm_1) ** self.m1
        # self.N_2 = 5 * 10 ** 6
        # self.sigm_2 = (self.N_D / self.N_2) ** (1 / self.m2) * self.sigm_d
        self.N_e = 10 ** 8
        self.sigm_e = (self.N_D / self.N_e) ** (1 / self.m2) * self.sigm_d

        plt.figure()
        x = [0, self.N_1, self.N_D, self.N_e]
        y = [self.sigm_1, self.sigm_1, self.sigm_d, self.sigm_e]
        plt.loglog(x, y, lw=2, marker='*')
        plt.xlabel('Cycle Numbers')
        plt.ylabel('Stress Amplitude/MPa')
        plt.xlim(10, 10 ** 8)
        plt.yticks([10, 100, 1000])
        plt.annotate(s='(%.1e,%.2f)' % (self.N_1, self.sigm_1), xy=(self.N_1, self.sigm_1))
        plt.annotate(s='(%.1e,%.2f)' % (self.N_D, self.sigm_d), xy=(self.N_D, self.sigm_d))
        plt.annotate(s='m1=%.2f' % self.m1, xy=(10 ** 3, 142))
        plt.annotate(s='m2=%.2f' % self.m2, xy=(10 ** 6.7, 40))
        plt.show()
        return


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Main()  # 创建QT对象
    window.show()  # QT对象显示
    sys.exit(app.exec_())
