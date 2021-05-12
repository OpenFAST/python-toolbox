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

