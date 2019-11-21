#!/usr/bin/env python
# _*_coding:utf-8 _*_
"""
@Time    : 10/30/19 1:35 PM
@Author  : Louis Wang
@Email   : louis1990.wang@outlook.com
@FileName: Fatigue_LouisW.py
@Software: PyCharm
"""

from numpy import empty, dot, hstack, delete, loadtxt
from xlrd import open_workbook
from collections import deque, defaultdict
from numba import jit
from numpy import array, roll, concatenate, argmax, abs


def loads_read(file_loc):
    """
	Reading the time series loads.
	The parameters of this function should modified according to real question now. This will fixed in the future.

	Args:
		file_loc: {str} File location of the load case occurrences numbers. xlsx file.

	Returns:
		loc_load: {dic} Dictionary with load case name: time series load
		loc_times: {dic} Dictionary with load case name: occurrence times.

	"""
    sheet = open_workbook(file_loc).sheet_by_index(0)
    loc_load = {}
    loc_times = {}
    # For fixed hun center load, take the [:,1:7] of the data from txt file.
    for loc_key, loc_name in enumerate(sheet.col_values(1)[1:]):
        loc_load[loc_name.split("\\")[-1]] = (
            loadtxt(r"%s.txt" % loc_name, skiprows=2)[:, 1:-2] * 1000
        )
        loc_times[loc_name.split("\\")[-1]] = (
            (sheet.col_values(7)[loc_key + 1] * 3600 / 600 * 20)
            if isinstance(sheet.col_values(7)[loc_key + 1], float)
            else sheet.col_values(3)[loc_key + 1]
        )

    return loc_load, loc_times


def unit_read(file_loc, load_num):
    """
		Read the unit load result and combined into a array. Only vensys file type supported now.

	Args:
		file_loc: Unit load file location.
		load_num: Number of the unit load result files

	Returns:
		Stress components of each node under unit load. /Nm-MPa

	"""
    ele_nd = deque()
    unit_load = {}
    for loadi in range(load_num):
        ele_s = deque()
        with open(r"%s/Sxyz_LF_%s.txt" % (file_loc, loadi + 1)) as f:
            line = f.readline()
            while line:
                try:
                    ele_s.append(
                        [
                            float(line[:8]),
                            float(line[8:22]),
                            float(line[22:35]),
                            float(line[35:]),
                        ]
                    )
                    line = f.readline()
                except ValueError:
                    line = f.readline()
            try:
                ele_nd = hstack((ele_nd, delete(array(ele_s), 0, 1)))
            except ValueError:
                ele_nd = array(ele_s)
    for loadi in range(len(ele_nd)):
        unit_load[ele_nd[loadi, 0]] = ele_nd[loadi, 1:]
    return unit_load


def _get_stress_function(d_method):
    """
	Stress combination function.

	Args:
		d_method: Method for stress calculation. max_principle, mises, max_shear, critical_plane

	Returns:
	Function for stress calculation.

    """
    if d_method == "max_principle":

        def func(x):
            return (x[:, 0] + x[:, 1]) / 2 + (
                (x[:, 0] - x[:, 1]) ** 2 / 4 + x[:, 2] ** 2
            ) ** (1 / 2)

    else:
        # Just return My
        def func(x):
            return x[:4]

    return func


def stress_combine(ts_load, u_load, d_method="max_principle"):
    """
	Stress combination according to the method used to calculate the fatigue damage.

	Args:
		ts_load: time series load of specific load case. KN-KNm
		u_load: unit load result of specific node. MPa-N
		d_method: the stress combination method used to calculate the fatigue damage.

	Returns:
		Time series stress data.

	"""
    func = _get_stress_function(d_method)
    stress_c = empty([len(ts_load), 3])
    #
    for i_load in range(3):
        stress_c[:, i_load] = dot(ts_load, u_load[0 + i_load:: 3])
    return func(stress_c)


def _get_peaks(series):
    """
    Return peaks & troughs of data, starting & finishing at the maximum

    Args:
        series: Initial time series data

    Returns:
        Spliced time series data started and end with the maximum value

    """
    # Eliminate repeated data points:
    series = series[series != roll(series, 1)]

    if len(series) == 0:
        return []

    # Split and rejoin at largest abs value
    max_idx = argmax(abs(series))
    series = concatenate((series[max_idx:], series[:max_idx]))

    # Find peaks and troughs
    prv = roll(series, -1)
    nxt = roll(series, 1)

    isPeak = ((series > prv) & (series > nxt)) | ((series < prv) & (series < nxt))

    # Close off the signal with the max (ie first) value, as required.
    return concatenate((series[isPeak], series[:1]))


def extract_cycles(series):
    """
    Iterate cycles in the series.

    Args:
        series: Iterable sequence of numbers

    Returns:
        tuple
		Each tuple contains three floats (low, high, times), where low and high
		define cycle amplitude and times equals to 1.0 for full cycles and 0.5
		for half cycles.

    """
    points = deque()
    for x in _get_peaks(series):
        points.append(x)
        while len(points) >= 3:
            # Form ranges X and Y from the three most recent points
            X = abs(points[-2] - points[-1])
            Y = abs(points[-3] - points[-2])

            if X < Y:
                # Read the next point
                break
            elif len(points) == 3:
                # Y contains the starting point
                # Count Y as one-half cycle and discard the first point
                yield points[0], points[1], 0.5
                points.popleft()
            else:
                # Count Y as one cycle and discard the peak and the valley of Y
                yield points[-3], points[-2], 1.0
                points[-3] = points[-1]
                points.pop()
                points.pop()
    else:
        # Count the remaining ranges as one-half cycles
        while len(points) > 1:
            yield points[0], points[1], 0.5
            points.popleft()


@jit(nopython=True)
def _round(x, ndigits=2):
    """Round cycle magnitudes to the given number of digits before counting"""
    return array([round(i, ndigits) for i in x])


def rainflow(series, ndigits=2):
    """
    Count cycles in the series.
    Args:
        series: Iterable sequence of numbers
        ndigits: Round cycle magnitudes to the given number of digits before counting

    Returns:
        A sorted list containing pairs of cycle magnitude, mean and count.
        One-half cycles are counted as 0.5, so the returned counts may not be whole numbers.

    """
    # Round cycle magnitudes to the given number of digits before counting
    series = _round(series, ndigits)
    counts = defaultdict(float)

    for low, high, times in extract_cycles(series):
        amplitude = abs(high - low) / 2
        mean = (high - low) / 2
        counts[(amplitude, mean)] += times
    return sorted(counts.items())


class FatigueFKM(object):
    """
	Fatigue damage calculation according to GL2010 and FKM.
	"""

    def __init__(
        self,
        mat="GGG",
        t=120,
        rm=360,
        rp=220,
        rz=125,
        r=-1,
        j_0=1,
        j=3,
        s_pu=2 / 3,
        gamma_m=1.265,
        modified_t=False,
    ):
        """
			Basic SN curve parameters calculation according to GL2010 Fig.5.B.3.
		Args:
			mat: {str} Material type. GGG for cast iron.
			rm: {float} Strength of tension. /MPa
			rp: {float} Strength of yield. /MPa
			t: {float} Thickness of the material.
			rz: {float} Surface roughness.
			r: {float} Stress ratio. default value is -1
			j: {float} Quality level for component.
			j_0: {float} Constant for material and test method.
			s_pu: {float} Survival probability.
			gamma_m: {float} Partial safety factor for material.
			modified_t: {bool} Whether to consider thickness correction. default is False. The thickness dependent
				material strength is considered in the tension/ yield strength
		Returns:
			Global parameters of the basic SN curve.
			M: {float} Mean stress sensitivity
			m1: {float} Slope of left of the knee point
			m2: {float} Slope of right of the knee point
			sigma_D: {float} Stress amplitude at the knee point
			N_D: {float} Cycle number of the knee point
		"""
        from math import log10, sqrt

        self.Rm = rm
        self.Rp = rp
        self.mat = mat
        self.t = t
        if modified_t:
            self.sign_tc = 1
        else:
            self.sign_tc = 0
        self.Rz = rz
        self.R = r
        self.j = j
        self.j_0 = j_0
        self.S_pu = s_pu
        self.gamma_M = gamma_m
        self.sigma_b = self.Rm * 1.06
        if self.mat == "GGG":
            self.M = 0.00035 * self.sigma_b + 0.08
        else:
            self.M = 0.00035 * self.sigma_b + 0.05
        # Surface roughness factor
        F_o = (
            1
            - 0.22 * (log10(self.Rz) ** 0.64) * log10(self.sigma_b)
            + 0.45 * (log10(self.Rz) ** 0.53)
        )
        # Notch factor and reduction factor
        alpha_k = 1
        n = 1
        belt_k = alpha_k / n
        F_ok = sqrt(belt_k ** 2 - 1 + 1 / (F_o ** 2))
        # Fatigue strength of specimen
        if self.mat == "GGG":
            sigma_w = 0.27 * self.sigma_b + 100
        else:
            sigma_w = 0.27 * self.sigma_b + 85
        # Fatigue strength of component
        sigma_wk = sigma_w / F_ok
        # Slopes of SN curve m1 and m2
        self.m1 = 5.5 / (F_ok ** 2) + 6
        self.m2 = 2 * self.m1 - 1
        # Factor for influence of mean stress. if R == -1: M = 1.

        if self.R == -1:
            Fm = 1
        else:
            a = (1 + self.R) / (1 - self.R) * sigma_wk / self.sigma_b
            u = 1 / (self.M + 1) * sigma_wk / self.sigma_b
            p = (1 / (self.M + 1) - 1 + u ** 2) / (u ** 2 - u)
            if p <= 1:
                Fm = -1 * (1 + p * a) / (2 * a ** 2 * (1 - p)) + sqrt(
                    1 / (1 - p) / a ** 2 + ((1 + p * a) / 2 / a ** 2 / (1 - p)) ** 2
                )
            else:
                Fm = -1 * (1 + p * a) / (2 * a ** 2 * (1 - p)) - sqrt(
                    1 / (1 - p) / a ** 2 + ((1 + p * a) / 2 / a ** 2 / (1 - p)) ** 2
                )
        # Stress amplitude at knee of SN curve
        sigma_kn = sigma_wk * Fm
        # Number of load cycles at knee of SN curve
        self.N_D = 10 ** (6.8 - 3.6 * (1 / self.m1))
        # Upgrading factors
        S_d = 0.85 ** (self.j - self.j_0)
        S_t = (self.t / 25) ** ((-0.15) * self.sign_tc)
        S = self.S_pu * S_d * S_t
        # Upgraded stress amplitude at knee of SN curve
        self.sigma_D = sigma_kn * S / self.gamma_M
        # Inflection point of the haigh diagram
        self.h1x = (
            self.sigma_D - self.M * self.sigma_D / (self.M - 1) - self.Rp / self.gamma_M
        )
        self.h2x = self.sigma_D / (self.M - 1)
        self.h3x = (self.Rp - self.sigma_D) / (1 - self.M)

        return

    def para_save(self, filename="SN_Curve"):
        """
		Save the parameters of basic sn curve and haigh diagram to .txt file.

		Args:
			filename: name of the file to be saved.

		Returns:
			filename.txt file
		"""
        filename = filename
        # Parameters of the haigh diagram, pulse stress amplitude.
        sigma_p = self.sigma_D / (self.M + 1)
        with open("%s.txt" % filename, "w") as f:
            f.write("Summary of the SN curve Parameter\n")
            f.write("*" * 30 + "\n")
            f.write("Material:\n%s\n" % self.mat)
            f.write("Stress ratio R:\n%s\n" % self.R)
            f.write("*" * 30 + "\n")
            f.write(
                "Calculation of SN curve according to the GL2010 5.B.3.1 and 5.B.3.2.\n"
            )
            f.write("*" * 30 + "\n")
            f.write("Tensile strength:\n")
            f.write("%s\n" % self.Rm)
            f.write("Yield Strength:\n")
            f.write("%s\n" % self.Rp)
            f.write("Alternating stress limit (Amplitude/Range):\n")
            f.write("%.2f/%.2f\n" % (self.sigma_D, 2 * self.sigma_D))
            f.write("Pulsating stress limit (Amplitude/Range):\n")
            f.write("%.2f/%.2f\n" % (sigma_p, 2 * sigma_p))
            f.write("Slope of the SN curve (Left/Right of the knee point):\n")
            f.write("%.2f/%.2f\n" % (self.m1, self.m2))
            f.write("Number of load cycles at knee of SN curve:\n")
            f.write("%.2e\n" % self.N_D)
            f.write("Mean stress sensitivity:\n")
            f.write("%.3f\n" % self.M)
            f.write("*" * 30 + "\n")
            f.write("Parameter of the Haigh Diagram:\n")
            f.write("*" * 30 + "\n")
            f.write("Mean stress / Amplitude:\n")
            f.write("%11.2f/%11d\n" % ((self.Rp / self.gamma_M), 0))
            if self.mat == "GGG":
                f.write(
                    "%11.2f/%11.2f\n"
                    % (
                        (self.Rp / self.gamma_M - self.sigma_D) / (1 - self.M),
                        self.M * (self.Rp / self.gamma_M - self.sigma_D) / (self.M - 1)
                        + self.sigma_D,
                    )
                )
            else:
                f.write("")
            f.write("%11.2f/%11.2f\n" % (sigma_p, sigma_p))
            f.write("%11d/%11.2f\n" % (0, self.sigma_D))
            f.write(
                "%11.2f/%11.2f\n"
                % (self.sigma_D / (self.M - 1), -self.sigma_D / (self.M - 1))
            )
            f.write("%11.2f/%11d\n" % ((-self.Rp / self.gamma_M), 0))
        return

    def sn_plot(self):
        """
			plot the basic SN curve
		"""
        import matplotlib.pyplot as plt

        # Plot of the basic SN curve according to GL2010
        sigma_1 = self.Rp * (1 - self.R) / self.gamma_M
        # Number of load cycles at upper fatigue limit
        N_1 = self.N_D * (2 * self.sigma_D / sigma_1) ** self.m1
        N_e = 10 ** 9
        sigma_e = (self.N_D / N_e) ** (1 / self.m2) * self.sigma_D
        x = [0, N_1, self.N_D, N_e]
        y = [sigma_1, sigma_1, self.sigma_D, sigma_e]
        plt.loglog(x, y, lw=2, marker="*")
        plt.xlabel("Cycle Numbers")
        plt.ylabel("Stress Amplitude/MPa")
        plt.xlim(10, 10 ** 9)
        plt.yticks([10, 100, 1000])
        plt.annotate(s="(%.2e,%.2f)" % (N_1, sigma_1), xy=(N_1, sigma_1))
        plt.annotate(
            s="(%.2e,%.2f)" % (self.N_D, self.sigma_D), xy=(self.N_D, self.sigma_D)
        )
        plt.annotate(s="m1=%.2f" % self.m1, xy=(10 ** 3, 142))
        plt.annotate(s="m2=%.2f" % self.m2, xy=(10 ** 7, 40))
        plt.show()
        return

    def factor_ms(self, mean_stress, ms_correction=True):
        """
			mean stress correction factor of the sn curve stress limit at knee point according to FKM.

		Args:
			mean_stress: {float} mean stress
			ms_correction: {Bool} mean stress correction.

		Returns:
			correction factor: sigma_d/ self.sigma_D
		"""
        if not ms_correction:
            return 1
        if self.h2x <= mean_stress <= self.h3x:
            return (self.sigma_D - self.M * mean_stress) / self.sigma_D
        elif self.h1x <= mean_stress < self.h2x:
            return (self.sigma_D - self.M * self.sigma_D / (self.M - 1)) / self.sigma_D
        elif mean_stress < self.h1x:
            return (mean_stress + self.Rp / self.gamma_M) / self.sigma_D
        else:
            return (
                self.Rp / self.gamma_M
                - (self.Rp / self.gamma_M - self.sigma_D) / (1 - self.M)
            ) / self.sigma_D

    def damage_node(self, markov, ms_correction=True):
        """
			Damage accumulation
		Args:
			markov: markov matrix. ([Amplitude, Mean}, Counts)
			ms_correction: Mean stress correction.

		Returns:
			Damage
		"""

        d = 0
        for i in iter(markov):
            amplitude = i[0][0]
            mean = i[0][1]
            count = i[1]
            if not amplitude == 0:
                sigma_d = self.sigma_D * self.factor_ms(mean, ms_correction)
                if mean <= sigma_d:
                    d += count / (self.N_D * (sigma_d / amplitude) ** self.m2)
                else:
                    d += count / (self.N_D * (sigma_d / amplitude) ** self.m1)
        return d

    def damage_cast(self, load_series, load_times, unit_load_result, node_list="all", ndigits=2):
        """
        Fatigue damage calculation for cast component.
        Args:
            load_series: Time series load, the format corresponds to the unit load stress result
            load_times: Times for time series load.
            unit_load_result: Stress results under unit load. This file is calculated
            node_list: Node list to be calculated of fatigue damage.
            ndigits: Round cycle magnitudes to the given number of digits before counting

        Returns:

        """
        # Dataframe to save the fatigue damage
        D_Details = DataFrame(
            empty([len(load_series), len(unit_load_result)]),
            columns=unit_load_result.keys(),
            index=load_series.keys(),
        )
        # List of the node and load case to be calculated.
        if node_list == "all":
            node_list = list(unit_load_result.keys())
        else:
            node_list = list(unit_load_result.keys())[node_list[0]: node_list[1]]
        loc_list = list(load_series.keys())
        d_max = [0, 0]
        # Loop for fatigue damage calculation
        for node, loc in product(node_list, loc_list):
            stress_h = stress_combine(load_series[loc], unit_load_result[node])
            D_Details.loc[loc, node] = self.damage_node(rainflow(stress_h, ndigits))*load_times[loc]
            if D_Details[node].sum() > d_max[1]:
                d_max = [node, D_Details[node].sum()]
            if node % 100 == 0:
                print(
                    "The maximum fatigue damage is %s for node %s"
                    % (d_max[0], d_max[1])
                )
        return D_Details


if __name__ == "__main__":
    import time
    from pandas import DataFrame
    from itertools import product

    # Fatigue damage calculation function
    test = FatigueFKM()

    # read the unit load result
    start = time.time()
    unit_result = unit_read(
        r"D:\Code\FatigueDamage\Project\3P5_Axle\Unit_Load_Results", 6
    )
    end = time.time()
    print("Reading unit load result finish. Time: %.2fs.\n" % (end - start))

    # read time series load
    start = time.time()
    load, times = loads_read(r"D:\Code\FatigueDamage\Data\times.xlsx")
    end = time.time()
    print("Reading time series load finish. Time: %.2fs.\n" % (end - start))

    # time series stress combination, rainflow and fatigue damage calculation.
    print("Start calculating of fatigue damage......\n")
    start = time.time()
    D = test.damage_cast(load, times, unit_result, [0, 100], 1)
    end = time.time()
    print("Fatigue damage calculation finish. Time: %.2fs.\n" % (end - start))
    D_Sum = D.apply(sum)
