"""
@Description:  SN curve calculation and modification according to GL2010.
@Author: Louis.Wang
@version: v1
@Date: 2019-09-23 13:55:31
@LastEditors: Louis.Wang
@LastEditTime: 2019-09-24 10:18:03
"""


class SN_curve():
    def __init__(self, mat=1, t=130, j=3, gamma_m=1.265, *args, **kwargs):
        # super().__init__(*args, **kwargs) 
        """
        @description:
            Basic paramater for SN curve calculation
        @para:
            mat: {int} Type of the material 1 for spherpodal graphire cast iron; 2 for cast steel.
            t: {float} Thickness of the material. /mm
        @return:
        """

        # Strength of tension. /MPa
        self.Rm = 360
        # strength of yield. /MPa
        self.Rp = 220
        # Type of the material 1 for spherpodal graphire cast iron; 2 for cast steel.
        self.mat = mat
        # Thickness of the material. /mm
        self.t = t
        # Thickness correction.  ///The thickness is not considered here
        self.sign_tc = 0
        # Surface roughness. 
        self.Rz = 125
        # stress ratio
        self.R = -1
        # quality level for component
        self.j = j
        # constant for material and test method  
        self.j_0 = 1
        # survival probability
        self.S_pu = 2 / 3
        # partial safety factor for material
        self.gamma_M = gamma_m
        # thickness correction. 0-No 1-Yes
        self.sign_t = 0

    def sn_base(self):
        """
        @description:
            Basic function for SN-curve paramater calculation according to the GL 2010.
        @param:
        @return:
            Parameters of the base SN curve.
        """
        from math import log10, sqrt

        self.Rb = self.Rm * 1.06
        self.Rs = self.Rp * 1.06

        # surface roughness factor
        F_o = 1 - 0.22 * (log10(self.Rz) ** 0.64) * log10(self.Rb) + 0.45 * (log10(self.Rz) ** 0.53)
        # notch factor
        alph_k = 1
        n = 1
        belt_k = alph_k / n
        # total influence factor fok
        F_ok = sqrt(belt_k ** 2 - 1 + 1 / (F_o ** 2))

        # fatigue strength of specimen
        if self.mat == 1:
            sigm_w = 0.27 * self.Rb + 100
            self.M = 0.00035 * self.Rb + 0.08
        else:
            sigm_w = 0.27 * self.Rb + 85
            self.M = 0.00035 * self.Rb + 0.05

        # fatigue shtrngth of component
        sigm_wk = sigm_w / F_ok

        # slopes of SN curve m1 and m2
        self.m1 = 5.5 / (F_ok ** 2) + 6
        self.m2 = 2 * self.m1 - 1

        # factor for influence of mean stress
        u = 1 / (self.M + 1) * sigm_wk / self.Rb
        a = (1 + self.R) / (1 - self.R) * sigm_wk / self.Rb
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
        # N_1*sigm_1**m = N_2*sigm_2**m
        self.N_1 = self.N_D * (2 * self.sigm_d / self.sigm_1) ** self.m1
        self.N_2 = 5 * 10 ** 6
        self.sigm_2 = (self.N_D / self.N_2) ** (1 / self.m2) * self.sigm_d
        self.N_e = 10 ** 9
        self.sigm_e = (self.N_2 / self.N_e) ** (1 / self.m2) * self.sigm_2

    def sn_base_plot(self):
        '''
        @description: 
            Figure plot of the initial SN surve.
        @param 
        @return: 
            Figure of the base SN curve.
        '''
        import matplotlib.pyplot as plt
        x = [0, self.N_1, self.N_D, self.N_2, self.N_e]
        y = [self.sigm_1, self.sigm_1, self.sigm_d, self.sigm_2, self.sigm_e]
        plt.loglog(x, y, lw=2, marker='*')
        plt.xlabel('Cycle Numbers')
        plt.ylabel('Stress Amplitude/MPa')
        plt.xlim(1, 10 ** 9)
        plt.yticks([10, 100, 1000])
        plt.annotate(s='(%.1e,%.2f)' % (self.N_1, self.sigm_1), xy=(self.N_1, self.sigm_1))
        plt.annotate(s='(%.1e,%.2f)' % (self.N_D, self.sigm_d), xy=(self.N_D, self.sigm_d))
        plt.annotate(s='(%.1e,%.2f)' % (5 * 10 ** 6, self.sigm_2), xy=(5 * 10 ** 6, 2 * self.sigm_2 - self.sigm_d))
        plt.annotate(s='m1=%.2f' % self.m1, xy=(10 ** 4, 150))
        plt.annotate(s='m2=%.2f' % self.m2, xy=(10 ** 7, 40))
        plt.show()
        return

    def haigh(self, plot=True):
        '''
        @description: 
            Parameter calculation of  Haigh digram.
        @param:
        @return: 
            Global Haigh diagram parameters.
        '''
        self.p1 = self.Rp / self.gamma_M  # 0
        if self.mat == 1:
            self.p2 = (self.Rp / self.gamma_M - self.sigm_d) / (
                        1 - self.M)  # (self.sigm_d-self.Rp*self.M/self.gamma_M)/(1-self.M)
        else:
            self.p2 = (self.Rp / self.gamma_M * 3 - self.sigm_d) / (
                        1 - self.M / 3)  # (self.sigm_d-self.Rp*self.M/3/self.gamma_M)/(1-self.M/3)
        self.p3 = self.sigm_d / (1 + self.M)  # self.sigm_d/(1+self.M)
        self.p4 = 0  # self.sigm_d
        self.p5 = self.sigm_d / (self.M - 1)  # self.sigm_d/(1-self.M)
        self.p6 = self.sigm_d / (1 - self.M) - self.Rp / self.gamma_M  # self.sigm_d/(1-self.M)
        self.p7 = -self.Rp / self.gamma_M  # 0
        return

    def Amp_Correct(self, mean_stress):
        '''
        @description: 
            sigm_d correction based on the Haigh digram.
        @param:
        @return: 
            Modified SN curve parameters: sigm_d.
        '''
        if mean_stress <= self.p7:
            sigm_d = 0
        elif mean_stress <= self.p6:
            sigm_d = mean_stress + self.Rp / self.gamma_M
        elif mean_stress <= self.p5:
            sigm_d = self.sigm_d / (1 - self.M)
        elif mean_stress <= self.p3:
            sigm_d = self.sigm_d - mean_stress * self.M
        elif mean_stress <= self.p2:
            if self.mat == 1:
                sigm_d = self.sigm_d - mean_stress * self.M
            elif self.mat == 2:
                sigm_d = -self.M * mean_stress / 3 + self.sigm_d * (3 + self.M) / 3 / (1 + self.M)
        elif mean_stress <= self.p1:
            sigm_d = self.Rp / self.gamma_M - mean_stress
        else:
            sigm_d = 0
        return sigm_d


if __name__ == "__main__":
    a = SN_curve(mat=1, t=130)
    a.sn_base()
    a.sn_base_plot()
