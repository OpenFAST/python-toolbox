import unittest
import os
import numpy as np
from pyFAST.converters.beam import ComputeStiffnessProps, ComputeInertiaProps, TransformCrossSectionMatrix


class Test(unittest.TestCase):

    def test_BD_H2(self):
        
        
        K = np.array([[ 3.22114734e+09, 3.25671541e+07, 0.00000000e+00, 0.00000000e+00,  0.00000000e+00, -4.99668423e+07],
              [ 3.25671541e+07, 2.43648638e+09, 0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  3.71784304e+07],
              [ 0.00000000e+00, 0.00000000e+00, 2.44325680e+10, 1.90703565e+08,  1.78086138e+09,  0.00000000e+00],
              [ 0.00000000e+00, 0.00000000e+00, 1.90703565e+08, 4.84065153e+10,  4.60061848e+09,  0.00000000e+00],
              [ 0.00000000e+00, 0.00000000e+00, 1.78086138e+09, 4.60061848e+09,  5.23688767e+10,  0.00000000e+00],
              [-4.99668423e+07, 3.71784304e+07, 0.00000000e+00, 0.00000000e+00,  0.00000000e+00,  2.33148898e+10]])

        I = np.array([[1439.88989372 ,    0.         ,    0.         ,    0.         ,    0.         ,   -5.92729626],
                    [   0.         , 1439.88989372 ,    0.         ,    0.         ,    0.         ,   -68.63720816],
                    [   0.         ,    0.         , 1439.88989372 ,    5.92729626 ,   68.63720816 ,    0.        ],
                    [   0.         ,    0.         ,    5.92729626 , 2594.38760678 ,   47.71380715 ,    0.        ],
                    [   0.         ,    0.         ,   68.63720816 ,   47.71380715 , 3434.04030538 ,    0.        ],
                    [  -5.92729626 ,  -68.63720816 ,    0.         ,    0.         ,    0.         ,   6028.42791216]])
        
        stiff = ComputeStiffnessProps()
        inertia = ComputeInertiaProps()
        transform = TransformCrossSectionMatrix()
        xs , ys = stiff.ComputeShearCenter(K)
        xt , yt = stiff.ComputeTensionCenter(K)
        xm , ym = inertia.ComputeMassCenter(I)

        # Approach BECAS
        # Translate to tension (elastic) center
        Kel = transform.CrossSectionRotoTranslationMatrix(K, xt, yt, 0.)
        # Find delta
        DeltaBecas = stiff.OrientationPrincipalAxesBecas(Kel)
        # Rotate by delta
        Kel_DeltaBecas = transform.CrossSectionRotoTranslationMatrix(Kel, 0., 0., DeltaBecas)

        # Approach ANBA4
        # Translate to tension (elastic) center
        Kdec = stiff.DecoupleStiffness(K)
        # Find delta
        DeltaANBA4 = stiff.PrincipalAxesRotationAngle(Kdec)
        # Rotate by delta
        Kel_DeltaANBA4 = transform.CrossSectionRotoTranslationMatrix(Kdec, 0., 0., -DeltaANBA4)

        print(stiff.ComputeTensionCenter(Kel_DeltaBecas))
        print(stiff.ComputeTensionCenter(Kel_DeltaANBA4))
        print(stiff.OrientationPrincipalAxesBecas(Kel_DeltaBecas))
        print(stiff.OrientationPrincipalAxesBecas(Kel_DeltaANBA4))
        print(stiff.PrincipalAxesRotationAngle(Kel_DeltaBecas))
        print(stiff.PrincipalAxesRotationAngle(Kel_DeltaANBA4))

        np.testing.assert_almost_equal(xs,  0.015468467117843322, 6)
        np.testing.assert_almost_equal(ys,  0.015668518364738194, 6)
        np.testing.assert_almost_equal(xt, -0.07288883346195946, 6)
        np.testing.assert_almost_equal(yt,  0.0078053017185913485, 6)
        np.testing.assert_almost_equal(xm, -0.04766837274110846, 6)
        np.testing.assert_almost_equal(ym,  0.004116492716458095, 6)
        np.testing.assert_almost_equal(DeltaBecas, -0.5874557755802033, 6)
        np.testing.assert_almost_equal(DeltaANBA4,  0.5874557755802033, 6)


if __name__ == '__main__':

    unittest.main()
