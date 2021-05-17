"""
Set of tools for stiffness and mass matrix of a beam section
"""

import numpy as np
from numpy import cos, sin

# --------------------------------------------------------------------------------}
# --- Bauchau 
# --------------------------------------------------------------------------------{
def K_sheartorsion_xbeam(J,K22,K23,K33,x2,x3):
    """ Returns Shear-torsion stiffness matrix.  See Eq.(13) of DyMore manual """
    return np.array( [
        [J + K22*x3**2 + 2*K23*x2*x3 + K33*x2**2, -K22*x3 - K23*x2, K23*x3 + K33*x2],
        [-K22*x3 - K23*x2, K22, -K23],
        [K23*x3 + K33*x2, -K23, K33]])

def K_axialbending_xbeam(S,H22,H23,H33,x2,x3):
    """ Returns Axial-Bending stiffness matrix. See Eq.(20) of DyMore manual """
    return np.array([
        [S, S*x3, -S*x2],
        [S*x3, H22 + S*x3**2, -H23 - S*x2*x3],
        [-S*x2, -H23 - S*x2*x3, H33 + S*x2**2]])

# --------------------------------------------------------------------------------}
# --- BeamDyn 
# --------------------------------------------------------------------------------{
def K_axialbending(EA, EI_x, EI_y, x_C=0, y_C=0, theta_p=0):
    """
    Axial bending problem. See KK for notations.
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    return np.array([
        [EA      , EA*y_C             , -EA*x_C            ] ,
        [EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C ] ,
        [-EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   ] 
        ])

def K_sheartorsion(GKt, GA, kxs, kys, x_S=0, y_S=0, theta_s=0):
    """
    Shear torsion problem. See KK for notations.
    """
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , -K_xx*y_S - K_xy*x_S                             ] ,
        [-K_xy                , K_yy                , K_xy*y_S + K_yy*x_S                              ] ,
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ]
        ])


def KK(EA, EI_x, EI_y, GKt, GA, kxs, kys, x_C=0, y_C=0, theta_p=0, x_S=0, y_S=0, theta_s=0):
    """ 
    Returns 6x6 stiffness matrix at the cross section origin O, based on inputs at centroid and shear center.
    INPUTS:
        - EA, EI_x, EI_y: diagonal terms for the axial bending expressed at the centroid and in the principal axis frame
        - GKt, GA*kxs, GA*kys: diagonal terms for the shear/torsion expressed at the shear center and in the princial shear direction frame
        - kxs, kys: dimensionless shear parameters
        - x_C, y_C: coordinates of the centroid (elastic center/ neutral axis), expressed FROM the origin of the cross section O
        - x_S, y_S:       "            shear center            "                  "                                             
        - theta_p : angle (around z) FROM the reference axes to the principal axes [rad]
        - theta_s :       "            "             "              principal shear axes [rad]
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , 0*EA    , 0*EA               , 0*EA               , -K_xx*y_S - K_xy*x_S                             ]    , 
        [-K_xy                , K_yy                , 0*EA    , 0*EA               , 0*EA               , K_xy*y_S + K_yy*x_S                              ]    , 
        [0*EA                 , 0*EA                , EA      , EA*y_C             , -EA*x_C            , 0*EA                                                ] , 
        [0*EA                 , 0*EA                , EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C , 0*EA                                                ] , 
        [0*EA                 , 0*EA                , -EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   , 0*EA                                                ] , 
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , 0*EA    , 0*EA               , 0*EA               , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ]
        ])

def MM(m,I_x,I_y,I_p,x_G=0,y_G=0,theta_i=0):
    """ 
    Returns the mass matrix at a given point O and with respect to given orientation axes based 
    on the values at the center of gravity and in the inertia axis frame.
    The convention is such that:
      - x_G,y_G      : the distaances FROM point O to point G
      - theta_i      : angle (around z) FROM the reference axes to the inertial axes
      - I_x, I_y, I_p: "diagonal" inertias for the body expressed in the inertial frame and at point G
    """
    Ixx = I_x*cos(theta_i)**2 + I_y*sin(theta_i)**2 
    Iyy = I_x*sin(theta_i)**2 + I_y*cos(theta_i)**2
    Ixy = (I_y-I_x)*sin(theta_i)*cos(theta_i)

    return np.array([
        [m      , 0*m     , 0*m      , 0*m                , 0*m                , -m*y_G]                    , 
        [0*m      , m     , 0*m      , 0*m                , 0*m                , m*x_G]                     , 
        [0*m      , 0*m     , m      , m*y_G            , -m*x_G           , 0*m]                         , 
        [0*m      , 0*m     , m*y_G  , Ixx + m*y_G**2   , -Ixy - m*x_G*y_G , 0*m]                         , 
        [0*m      , 0*m     , -m*x_G , -Ixy - m*x_G*y_G , Iyy + m*x_G**2   , 0*m]                         , 
        [-m*y_G , m*x_G , 0*m      , 0*m                , 0*m                , I_p + m*x_G**2 + m*y_G**2]
        ])

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
        M2 = T.T @ M1 @ T 
        # Rotation 
        R = self.CrossSectionRotationMatrix(alpha)
        M3 = R @ M2 @ R.T
        return M3

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

