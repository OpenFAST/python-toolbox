import numpy as np

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

class TransformCrossSectionMatrix(object):

    def CrossSectionTranslationMatrix(self, x, y):

        T = np.eye(6)
        T[0,5] = y
        T[1,5] = -x
        T[2,3] = -y
        T[2,4] = x

        return T

    def CrossSectionRotationMatrix(self, alpha):

        c=np.cos(alpha)
        s=np.sin(alpha)

        R1=[[c,s,0],
            [-s,c,0],
            [0,0,1]]
        R=np.vstack((np.hstack((R1, np.zeros((3,3)))),
           np.hstack((np.zeros((3,3)), R1))))

        return R

    def CrossSectionRotoTranslationMatrix(self, M1, x, y, alpha):
        
        # Translation
        T = self.CrossSectionTranslationMatrix(x, y)
        M2 = np.matmul(T.T, M1)
        M3 = np.matmul(M2, T) 

        # Rotation 
        R = self.CrossSectionRotationMatrix(alpha)
        M4 = np.matmul(np.matmul(R, M3), R.T)

        return M4
        
    def trsf_sixbysix(M, T):
        """
        Transform six-by-six compliance/stiffness matrix. 
        change of reference frame in engineering (or Voigt) notation.
        
        Parameters
        ----------
        M : np.ndarray
            6x6 Siffness or Mass Matrix
        T : np.ndarray
            Transformation Matrix
            
        Returns
        ----------
        res : np.ndarray
            Transformed 6x6 matrix
        """

        TS_1 = np.dot(np.dot(T.T, M[0:3, 0:3]), T)
        TS_2 = np.dot(np.dot(T.T, M[3:6, 0:3]), T)
        TS_3 = np.dot(np.dot(T.T, M[0:3, 3:6]), T)
        TS_4 = np.dot(np.dot(T.T, M[3:6, 3:6]), T)

        tmp_1 = np.vstack((TS_1, TS_2))
        tmp_2 = np.vstack((TS_3, TS_4))
        res = np.hstack((tmp_1, tmp_2))
        return res

class ComputeStiffnessProps(object):

    def ComputeShearCenter(self, K):   # shear center equiv. to elastic axes
        K1 = np.array([[K[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[K[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(K1, -K3)
        return [-Y[1,2], Y[0,2]]

    def ComputeTensionCenter(self, K):  # tension center equiv. to neutral axes
        K1 = np.array([[K[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[K[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(K1, -K3)
        return [Y[2,1], -Y[2,0]]

    def OrientationPrincipalAxesBecas(self, K):
        
        ksub=K[3:5,3:5]
        [ val, mod ] = np.linalg.eig(ksub)
        val = np.sort(np.diag(val))
        ind = np.argsort(np.diag(val))
        mod = mod[:,ind]
        Delta = np.arctan(mod[1,0]/mod[0,0])
        
        return Delta

    def DecoupleStiffness(self, K):
        K1 = np.array([[K[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[K[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(K1, -K3)
        I3 = np.eye(3)
        Z3 = np.zeros((3,3))
        TL = np.block([[I3, Z3], [Y.T, I3]])
        TR = np.block([[I3, Y], [Z3, I3]])
        return TL @ K @ TR

    def PrincipalAxesRotationAngle(self, decoupledK):
        K1 = np.array([[decoupledK[i, j] for j in range(3)] for i in range(3)])
        K3 = np.array([[decoupledK[i+3, j+3] for j in range(3)] for i in range(3)])
        (w1, v1) = np.linalg.eig(K1)
        (w3, v3) = np.linalg.eig(K3)
        
        if np.abs(v3[0,0]) < np.abs(v3[0,1]):
            angle = np.arccos(v3[0,0])
        else:
            angle = -np.arcsin(v3[0,1])
        return angle

class ComputeInertiaProps(object):     
    
    def ComputeMassCenter(self, M):
        M1 = np.array([[M[i, j] for j in range(3)] for i in range(3)])
        M3 = np.array([[M[i, j+3] for j in range(3)] for i in range(3)])
        Y = np.linalg.solve(M1, -M3)
        return [Y[2,1], -Y[2,0]]

if __name__ == "__main__":

    stiff = ComputeStiffnessProps()
    inertia = ComputeInertiaProps()
    transform = TransformCrossSectionMatrix()
    xe , ye = stiff.ComputeShearCenter(K)
    xt , yt = stiff.ComputeTensionCenter(K)
    xm , ym = inertia.ComputeMassCenter(I)

    # Approach BECAS
    Kel = transform.CrossSectionRotoTranslationMatrix(K, xe, ye, 0.)
    DeltaBecas = stiff.OrientationPrincipalAxesBecas(Kel)

    # Approach ANBA4
    Kdec = stiff.DecoupleStiffness(K)
    DeltaANBA4 = stiff.PrincipalAxesRotationAngle(Kdec)

if __name__ == '__test__':

    stiff = ComputeStiffnessProps()
    inertia = ComputeInertiaProps()
    transform = TransformCrossSectionMatrix()
    xe , ye = stiff.ComputeShearCenter(K)
    xt , yt = stiff.ComputeTensionCenter(K)
    xm , ym = inertia.ComputeMassCenter(I)

    # Approach BECAS
    Kel = transform.CrossSectionRotoTranslationMatrix(K, xe, ye, 0.)
    DeltaBecas = stiff.OrientationPrincipalAxesBecas(Kel)

    # Approach ANBA4
    Kdec = stiff.DecoupleStiffness(K)
    DeltaANBA4 = stiff.PrincipalAxesRotationAngle(Kdec)

    np.testing.assert_almost_equal(xe,  0.015468467117843322, 6)
    np.testing.assert_almost_equal(ye,  0.015668518364738194, 6)
    np.testing.assert_almost_equal(xt, -0.07288883346195946, 6)
    np.testing.assert_almost_equal(yt,  0.0078053017185913485, 6)
    np.testing.assert_almost_equal(xm, -0.04766837274110846, 6)
    np.testing.assert_almost_equal(ym,  0.004116492716458095, 6)
    np.testing.assert_almost_equal(DeltaBecas, -0.5780573896723843, 6)
    np.testing.assert_almost_equal(DeltaANBA4,  0.5874557755802033, 6)


