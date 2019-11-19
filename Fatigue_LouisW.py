#!/usr/bin/env python
# _*_coding:utf-8 _*_
"""
@Time    : 10/30/19 1:35 PM
@Author  : Louis Wang
@Email   : louis1990.wang@outlook.com
@FileName: Fatigue_LouisW.py
@Software: PyCharm
"""

from numpy import empty, dot, array, hstack, delete, loadtxt, around, concatenate
from scipy.signal import butter, filtfilt
from xlrd import open_workbook
from collections import deque, defaultdict
from numba import jit
import functools
import threading


def loads_read(file_loc):
    """
    Reading the time series loads.
    This function in only support for the vensys load file type now. Support for goldwind load file type will added in
    the future.

    Args:
        file_loc: {str} File location of the load case occurrences numbers. xlsx file.

    Returns:
        loc_load: {dic} Dictionary with load case name: time series load
        loc_times: {dic} Dictionary with load case name: occurrence times.

    """
    sheet = open_workbook(file_loc).sheet_by_index(0)
    loc_load = {}
    loc_times = {}

    for loc_key, loc_name in enumerate(sheet.col_values(1)[1:]):
        loc_load[loc_name.split("\\")[-1]] = loadtxt(r"%s.txt" % loc_name, skiprows=2)[
            :, -6:-1
        ]
        loc_times[loc_name.split("\\")[-1]] = (
            (sheet.col_values(7)[loc_key + 1] * 3600 / 600 * 20)
            if isinstance(sheet.col_values(7)[loc_key + 1], float)
            else sheet.col_values(3)[loc_key + 1]
        )
        ########################################################
        # Expend the blade root load to three blade roots load.
        # This part should modified according to actual conditions.
        loc_load[loc_name.split("\\")[-1]] = (
            concatenate(
                (
                    loc_load[loc_name.split("\\")[-1]],
                    loc_load[loc_name.split("\\")[-1]],
                    loc_load[loc_name.split("\\")[-1]],
                ),
                axis=1,
            )
            * 1000
        )
        ########################################################
    return loc_load, loc_times


def unit_read(file_loc, load_num):
    """
        Read the unit load result and combined into a array. Only vensys file type supported now.

    Args:
        file_loc: unit load file location.
        load_num: number of the unit load result files

    Returns:
        ele_nd: {array} node number and stress component under unit load. /Nm-MPa

    """
    ele_nd = []
    unit_load = {}
    for loadi in range(load_num):
        ele_s = []
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


def reversals(series, left=False, right=False):
    """
    Iterate reversal points in the series.
    A reversal point is a point in the series at which the first derivative
    changes sign. Reversal is undefined at the first (last) point because the
    derivative before (after) this point is undefined. The first and the last
    points may be treated as reversals by setting the optional parameters
    `left` and `right` to True.
    Parameters
    ----------
    series : iterable sequence of numbers
    left: bool, optional
        If True, yield the first point in the series (treat it as a reversal).
    right: bool, optional
        If True, yield the last point in the series (treat it as a reversal).
    Yields
    ------
    float
        Reversal points.
    """
    series = iter(series)

    x_last, x = next(series), next(series)
    d_last = x - x_last

    if left:
        yield x_last
    for x_next in series:
        if x_next == x:
            continue
        d_next = x_next - x
        if d_last * d_next < 0:
            yield x
        x_last, x = x, x_next
        d_last = d_next
    if right:
        yield x_next


def _sort_lows_and_highs(func):
    """Decorator for extract_cycles"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        for low, high, times in func(*args, **kwargs):
            if low < high:
                yield low, high, times
            else:
                yield high, low, times

    return wrapper


@_sort_lows_and_highs
def extract_cycles(series, left=False, right=False):
    """Iterate cycles in the series.
    Parameters
    ----------
    series : iterable sequence of numbers
    left: bool, optional
        If True, treat the first point in the series as a reversal.
    right: bool, optional
        If True, treat the last point in the series as a reversal.
    Yields
    ------
    cycle : tuple
        Each tuple contains three floats (low, high, times), where low and high
        define cycle amplitude and times equals to 1.0 for full cycles and 0.5
        for half cycles.
    """
    points = deque()

    for x in reversals(series, left=left, right=right):
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
                last = points.pop()
                points.pop()
                points.pop()
                points.append(last)
    else:
        # Count the remaining ranges as one-half cycles
        while len(points) > 1:
            yield points[0], points[1], 0.5
            points.popleft()


def rainflow(series, left=True, right=True):
    """Count cycles in the series.
        Parameters
        ----------
            series : iterable sequence of numbers
            left: bool, optional
                If True, treat the first point in the series as a reversal.
            right: bool, optional
                If True, treat the last point in the series as a reversal.
        Returns
        -------
            A sorted list containing pairs of cycle magnitude, mean and count.
            One-half cycles are counted as 0.5, so the returned counts may not be
            whole numbers.
        """
    counts = defaultdict(float)
    for low, high, times in extract_cycles(series, left=left, right=right):
        delta = abs(high - low) / 2
        mean = (high - low) / 2
        counts[(delta, mean)] += times
    return sorted(counts.items())


def stress_combine(
    ts_load,
    u_load,
    d_method="max_principle",
    filtering=True,
    ndigits=2,
    order=8,
    cutoff=0.3,
):
    """
        Stress combination according to the method used to calculate the fatigue damage.
        Digital filtering of the time series stress history data (low pass).
        Signal frequency fs<2Hz fow wind load.
        Sampling frequency fm=20Hz for vensys, 50Hz for goldwind.
        Cut off frequency fc=3Hz.

    Args:
        ts_load: time series load. KN-KNm
        u_load: unit load result. MPa-N
        d_method: the method used to calculate the fatigue damage.
        filtering: {bool} Load signal filtering. Detail parameters described in signal_filter function.
        ndigits: {int} Number of decimal places to round to (default: 2)
        order: {int} Order of signal filter.
        cutoff: {float}  By default is 0.3 = 2*fc/fm for vensys.

    Returns:
        Filtered and precision controlled time series stress history data.

    """
    stress_c = {}
    stress_r = {}

    for key_load in ts_load:
        loc_num = int(len(ts_load[key_load]))
        stress_c[key_load] = empty([loc_num, 3])
        for i_load in range(3):
            stress_c[key_load][:, i_load] = dot(ts_load[key_load], u_load[0 + i_load:: 3])

    if d_method == "max_principle":
        stress_r[key_load] = (
            stress_c[key_load][:, 0] + stress_c[key_load][:, 1]
        ) / 2 + (
            (stress_c[key_load][:, 0] - stress_c[key_load][:, 1]) ** 2 / 4
            + stress_c[key_load][:, 2] ** 2
        ) ** (
            1 / 2
        )
    else:
        # More function/fatigue damage accumulation hypothesis will be added in the future
        stress_r[key_load] = stress_c[key_load][:, 4].reshape([loc_num, 1])
    b, a = butter(order, cutoff)
    return (
        [
            around(filtfilt(b, a, stress_r[key_stress]), ndigits)
            for key_stress in stress_r
        ]
        if filtering
        else [around(stress_r[key_stress], ndigits) for key_stress in stress_r]
    )


def chunks(arr, m):
    """
    Divide a list to m parts.

    Args:
        arr: List to be divided
        m: Number of parts (Number of multi-processing).

    Returns:
        {list} Divided list data.

    """
    n = round(len(arr) / (m - 1))
    return [arr[id_node : id_node + n] for id_node in range(0, len(arr), n)]


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
            sigma_d_m = self.sigma_D - self.M * mean_stress
        elif self.h1x <= mean_stress < self.h2x:
            sigma_d_m = self.sigma_D - self.M * self.sigma_D / (self.M - 1)
        elif mean_stress < self.h1x:
            sigma_d_m = mean_stress + self.Rp / self.gamma_M
        else:
            sigma_d_m = self.Rp / self.gamma_M - (
                self.Rp / self.gamma_M - self.sigma_D
            ) / (1 - self.M)
        return sigma_d_m / self.sigma_D

    def damage_s(self, amplitude, mean, count, ms_correction=True):
        """
			Damage accumulation
		Args:
			amplitude:
			mean:
			count:
			ms_correction:

		Returns:
			Damage
		"""
        if amplitude == 0:
            return 0
        sigma_d = self.sigma_D * self.factor_ms(mean, ms_correction)
        if mean <= sigma_d:
            damage = count / (self.N_D * (sigma_d / amplitude) ** self.m2)
        else:
            damage = count / (self.N_D * (sigma_d / amplitude) ** self.m1)
        return damage


class MyThread(threading.Thread):

    """
    Multi-Processing
    """

    def __init__(self, mid, nodelist):
        threading.Thread.__init__(self)
        self.mid = mid
        self.node_list = nodelist

    def run(self):
        # 把要执行的代码写到run函数里面 线程在创建后会直接运行run函数
        print("Starting Threading-%s\n" % self.mid)

        # count_node = 1
        for node in self.node_list:
            stress_h = {}
            stress_h = stress_combine(
                load_series, unit_load_result[node], filtering=True, ndigits=2
            )
            #     for key in stress_h:
            #         markov = {}
            #         markov = rainflow(stress_h[key].ravel())
            #         for i in markov:
            #             D_Details.loc[key, node] += (
            #                 test.damage_s(i[0][0], i[0][1], i[1]) * load_times[key]
            #             )
            #
            # if count_node % 10 == 0:
            #     print(
            #         "%s/%s node damage calculation finish. Time: %.2fs."
            #         % (count_node, len(unit_load_result), (end - start))
            #     )
            # count_node += 1
        # D_Sum = D_Details.apply(sum)
        print("Exiting: Threading-%s\n" % self.mid)


if __name__ == "__main__":
    import time
    from pandas import DataFrame

    # Fatigue damage calculation function
    test = FatigueFKM()
    mul_threading = 8

    # read the unit load result
    start = time.time()
    unit_load_result = unit_read(r"D:\Code\FatigueDamage\Factor", 15)
    end = time.time()
    print("Reading unit load result finish. Time: %.2fs." % (end - start))

    # read time series load
    start = time.time()
    load_series, load_times = loads_read(r"D:\Code\FatigueDamage\Data\times.xlsx")
    end = time.time()
    print("Reading time series load finish. Time: %.2fs." % (end - start))

    # time series stress combination, rainflow and fatigue damage calculation.
    print(
        "Start calculating for fatigue damage with %s-threading process......"
        % mul_threading
    )
    node_list = chunks(list(unit_load_result.keys())[:100], mul_threading)
    start = time.time()
    # Dataframe to save the fatigue damage
    D_Details = DataFrame(
        empty([len(load_series), len(unit_load_result)]),
        columns=unit_load_result.keys(),
        index=load_series.keys(),
    )

    threads = {}
    for i in range(len(node_list)):
        threads[i] = MyThread(i, node_list[i])
    for i in range(len(node_list)):
        threads[i].start()
    for i in range(len(node_list)):
        threads[i].join()
    end = time.time()
    print("Fatigue damage calculation finish. Time: %.2fs." % (end - start))
