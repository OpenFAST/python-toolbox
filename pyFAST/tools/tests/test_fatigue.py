import unittest
import numpy as np
from pyFAST.tools.fatigue import *

class TestFatigue(unittest.TestCase):

    def test_leq_1hz(self):
        """Simple test of wetb.fatigue.eq_load using a sine
        signal.
        """
        amplitude = 1
        m = 1
        point_per_deg = 100

        for amplitude in [1,2,3]:
            peak2peak = amplitude * 2
            # sine signal with 10 periods (20 peaks)
            nr_periods = 10
            time = np.linspace(0, nr_periods*2*np.pi, point_per_deg*180)
            neq = time[-1]
            # mean value of the signal shouldn't matter
            signal = amplitude * np.sin(time) + 5
            r_eq_1hz = eq_load(signal, no_bins=1, m=m, neq=neq)[0]
            r_eq_1hz_expected = ((2*nr_periods*amplitude**m)/neq)**(1/m)
            np.testing.assert_allclose(r_eq_1hz, r_eq_1hz_expected)

            # sine signal with 20 periods (40 peaks)
            nr_periods = 20
            time = np.linspace(0, nr_periods*2*np.pi, point_per_deg*180)
            neq = time[-1]
            # mean value of the signal shouldn't matter
            signal = amplitude * np.sin(time) + 9
            r_eq_1hz2 = eq_load(signal, no_bins=1, m=m, neq=neq)[0]
            r_eq_1hz_expected2 = ((2*nr_periods*amplitude**m)/neq)**(1/m)
            np.testing.assert_allclose(r_eq_1hz2, r_eq_1hz_expected2)

            # 1hz equivalent should be independent of the length of the signal
            np.testing.assert_allclose(r_eq_1hz, r_eq_1hz2)

    def test_rainflow_combi(self):
        # Signal with two frequencies and amplitudes
        amplitude = 1
        # peak2peak = amplitude * 2
        m = 1
        point_per_deg = 100

        nr_periods = 10
        time = np.linspace(0, nr_periods*2*np.pi, point_per_deg*180)

        signal = (amplitude*np.sin(time)) + 5 + (amplitude*0.2*np.cos(5*time))
        cycles, ampl_bin_mean, ampl_edges, mean_bin_mean, mean_edges = \
            cycle_matrix(signal, ampl_bins=10, mean_bins=5)

        cycles.sum()



    def test_astm1(self):

        signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])

        ampl, mean = rainflow_astm(signal)
        np.testing.assert_array_equal(np.histogram2d(ampl, mean, [6, 4])[0], np.array([[ 0., 1., 0., 0.],
                                                                                                           [ 1., 0., 0., 2.],
                                                                                                           [ 0., 0., 0., 0.],
                                                                                                           [ 0., 0., 0., 1.],
                                                                                                           [ 0., 0., 0., 0.],
                                                                                                           [ 0., 0., 1., 2.]]))

    def test_windap1(self):
        signal = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
        ampl, mean = rainflow_windap(signal, 18, 2)
        np.testing.assert_array_equal(np.histogram2d(ampl, mean, [6, 4])[0], np.array([[ 0., 0., 1., 0.],
                                                                                       [ 1., 0., 0., 2.],
                                                                                       [ 0., 0., 0., 0.],
                                                                                       [ 0., 0., 0., 1.],
                                                                                       [ 0., 0., 0., 0.],
                                                                                       [ 0., 0., 2., 1.]]))

    def test_eq_load_basic(self):
        import numpy.testing
        signal1 = np.array([-2.0, 0.0, 1.0, 0.0, -3.0, 0.0, 5.0, 0.0, -1.0, 0.0, 3.0, 0.0, -4.0, 0.0, 4.0, 0.0, -2.0])
        try:
            M1=eq_load(signal1, no_bins=50, neq=[1, 17], m=[3, 4, 6], rainflow_func=rainflow_windap)
            doTest=True
        except FloatingPointError as e:
            doTest=False
            print('>>> Floating point error')
        M1_ref=np.array([[10.348414123746581, 9.635653414943068, 9.122399471334054], [4.024613313976801, 4.745357541147315, 5.68897815218057]])
        #M1_ref=np.array([[10.311095426959747, 9.5942535021382174, 9.0789213365013932],[4.010099657859783, 4.7249689509841746, 5.6618639965313005]])
        numpy.testing.assert_almost_equal(M1,M1_ref,decimal=5)
        #signal2 = signal1 * 1.1
        #         print (eq_load(signal1, no_bins=50, neq=17, rainflow_func=rainflow_windap))
        #         print (eq_load(signal1, no_bins=50, neq=17, rainflow_func=rainflow_astm))
        #         # equivalent load for default wohler slopes
        #         # Cycle matrix with 4 amplitude bins and 4 mean value bins
        #         print (cycle_matrix(signal1, 4, 4, rainflow_func=rainflow_windap))
        #         print (cycle_matrix(signal1, 4, 4, rainflow_func=rainflow_astm))
        #         # Cycle matrix where signal1 and signal2 contributes with 50% each
        #         print (cycle_matrix([(.5, signal1), (.5, signal2)], 4, 8, rainflow_func=rainflow_astm))



if __name__ == '__main__':
    unittest.main()
