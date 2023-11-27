import unittest
import numpy as np
import os
from openfast_toolbox.aeroacoustics import *

scriptDir = os.path.dirname(__file__)
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
class TestAeroAcousticTools(unittest.TestCase):

    def test_dummy(self):
        x=0
        np.testing.assert_equal(x, 0)

if __name__ == '__main__':
    unittest.main()
