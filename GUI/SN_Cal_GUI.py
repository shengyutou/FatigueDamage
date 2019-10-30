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
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtWidgets import QMessageBox
from math import log10, sqrt
import SNCurve

MainWindow_main = SNCurve.Ui_SN_Curve


class Main(QtWidgets.QMainWindow, MainWindow_main):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        MainWindow_main.__init__(self)
        self.setupUi(self)
        self.setWindowIcon(QtGui.QIcon("SN.png"))
        self.Cal_Plot.clicked.connect(self.sn_cal)

    def sn_cal(self):
        """
        Calculation of the SN curve
        :return:
        """
        # Para Read
        try:
            mat = str(self.l_mat.currentText())
            rp = float(self.l_Rp.text())
            rm = float(self.l_Rm.text())
            t = float(self.l_t.text())
            rz = float(self.l_Rz.text())
            j = float(self.l_j.text())
            j_0 = float(self.l_j_0.text())
            S_pu = float(self.l_S_pu.text())
            gamma_m = float(self.l_gamma_M.text())
            sigma_b = rm * 1.06
            R = -1
            if self.check_t_c.isChecked():
                sigma_t = 1
            else:
                sigma_t = 0
        except Exception as e:
            with open("ErrorLog.txt", "a") as f:
                f.write(
                    time.strftime("\n%Y-%m-%d-%H-%M", time.localtime(time.time()))
                    + "\n"
                )
                f.write(repr(e))
                f.close()
                QMessageBox.warning(self, "Error", "参数输入错误，请检查！")
            return

        # Calculation
        if mat == "GGG":
            M = 0.00035 * sigma_b + 0.08
        else:
            M = 0.00035 * sigma_b + 0.05

        # surface roughness factor
        F_o = (
            1 - 0.22 * (log10(rz) ** 0.64) * log10(sigma_b) + 0.45 * (log10(rz) ** 0.53)
        )

        # notch factor
        alpha_k = 1
        n = 1
        belt_k = alpha_k / n

        # total influence factor fok
        F_ok = sqrt(belt_k ** 2 - 1 + 1 / (F_o ** 2))

        # fatigue strength of specimen
        if mat == "GGG":
            sigma_w = (0.27 * sigma_b) + 100
        else:
            sigma_w = 0.27 * sigma_b + 85

        # fatigue strength of component
        sigma_wk = sigma_w / F_ok

        # slopes of SN curve m1 and m2
        m1 = 5.5 / (F_ok ** 2) + 6
        m2 = 2 * m1 - 1

        # factor for influence of mean stress
        u = 1 / (M + 1) * sigma_wk / sigma_b
        a = (1 + R) / (1 - R) * sigma_wk / sigma_b
        p = (1 / (M + 1) - 1 + u ** 2) / (u ** 2 - u)
        if a == 0:
            Fm = 1
        else:
            if p <= 1:
                Fm = -1 * (1 + p * a) / (2 * a ** 2 * (1 - p)) + sqrt(
                    1 / (1 - p) / a ** 2 + ((1 + p * a) / 2 / a ** 2 / (1 - p)) ** 2
                )
            else:
                Fm = -1 * (1 + p * a) / (2 * a ** 2 * (1 - p)) - sqrt(
                    1 / (1 - p) / a ** 2 + ((1 + p * a) / 2 / a ** 2 / (1 - p)) ** 2
                )

        # stress amplitude at knee of SN curve
        sigma_A = sigma_wk * Fm

        # number of load cycles at knee of SN curve
        N_D = 10 ** (6.8 - 3.6 * (1 / m1))

        # upgrading factors
        # quality level
        S_d = 0.85 ** (j - j_0)
        # thickness-dependent tensile value rm
        S_t = (t / 25) ** ((-0.15) * sigma_t)
        # total upgrading factor
        S = S_pu * S_d * S_t

        # upgraded stress amplitude at knee of SN curve
        sigma_d = sigma_A * S / gamma_m
        sigma_p = sigma_d / (M + 1)
        # upper limit of fatigue life line
        sigma_1 = rp * (1 - R) / gamma_m
        # number of load cycles at upper fatigue limit
        N_1 = N_D * (2 * sigma_d / sigma_1) ** m1
        N_e = 10 ** 8
        sigma_e = (N_D / N_e) ** (1 / m2) * sigma_d

        # SN curve parameter written to txt file.
        with open("SN_Curve.txt", "w") as f:
            f.write("Summary of the SN curve Parameter\n")
            f.write("Material:" + mat + "\n")
            f.write("Stress ratio R:" + str(R) + "\n")
            f.write("*" * 30 + "\n")
            f.write(
                "Calculation of the SN curve according to the GL2010 5.B.3.1 and 5.B.3.2.\n"
            )
            f.write("Tensile strength:\n")
            f.write(str(rm) + "\n")
            f.write("Yield Strength:\n")
            f.write(str(rp) + "\n")
            f.write("Alternating stress limit (Amplitude/Range):\n")
            f.write(str(round(sigma_d, 3)) + "/" + str(round(2 * sigma_d, 3)) + "\n")
            f.write("Pulsating stress limit (Amplitude/Range):\n")
            f.write(str(round(sigma_p, 3)) + "/" + str(round(2 * sigma_p, 3)) + "\n")
            f.write("Slope of the SN curve (Left/Right of the knee point):\n")
            f.write(str(round(m1, 3)) + "/" + str(round(m2, 3)) + "\n")
            f.write("Number of load cycles at knee of SN curve:\n")
            f.write(str("%.3e" % N_D) + "\n")
            f.write("Mean stress sensitivity:\n")
            f.write(str(round(M, 3)) + "\n")
            f.write("*" * 30 + "\n")
            f.write("Parameter of the Haigh Diagram:\n")
            f.write("Mean stress / Amplitude:\n")
            f.write(str(round(rp / gamma_m, 3)) + "/" + "0\n")
            if mat == "GGG":
                f.write(
                    str(round((rp / gamma_m - sigma_d) / (1 - M), 3))
                    + "/"
                    + str(round(M * (rp / gamma_m - sigma_d) / (M - 1) + sigma_d, 3))
                    + "\n"
                )
            else:
                f.write("")
            f.write(str(round(sigma_p, 3)) + "/" + str(round(sigma_p, 3)) + "\n")
            f.write("0" + "/" + str(round(sigma_d, 3)) + "\n")
            f.write(
                str(round(sigma_d / (M - 1), 3))
                + "/"
                + str(round(-sigma_d / (M - 1), 3))
                + "\n"
            )
            f.write(str(round(-rp / gamma_m, 3)) + "/" + "0\n")
            f.close()
        # Plot of the SN curve
        x = [0, N_1, N_D, N_e]
        y = [sigma_1 / 2, sigma_1 / 2, sigma_d, sigma_e]
        plt.loglog(x, y, lw=2, marker="*")
        plt.xlabel("Cycle Numbers")
        plt.ylabel("Stress Amplitude/MPa")
        plt.xlim(10, 10 ** 8)
        plt.yticks([10, 100, 500])
        plt.annotate(s="(%.1e,%.2f)" % (N_1, sigma_1 / 2), xy=(N_1, sigma_1 / 2))
        plt.annotate(s="(%.1e,%.2f)" % (N_D, sigma_d), xy=(N_D, sigma_d))
        plt.annotate(s="m1=%.2f" % m1, xy=(10 ** 3, 142))
        plt.annotate(s="m2=%.2f" % m2, xy=(10 ** 6.7, 40))
        plt.show()
        return


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = Main()  # 创建QT对象
    window.show()  # QT对象显示
    sys.exit(app.exec_())
