import numpy as np
from numpy import cos, sin
import pandas as pd
import os
# from weio.hawc2_htc_file import HAWC2HTCFile
# from weio.csv_file import CSVFile
# from weio.fast_input_file import FASTInputFile
from pyFAST.input_output.hawc2_htc_file import HAWC2HTCFile
from pyFAST.input_output.csv_file import CSVFile
from pyFAST.input_output.fast_input_file import FASTInputFile

from .beam import *

# --------------------------------------------------------------------------------}
# ---beamDynToHawc2
# --------------------------------------------------------------------------------{
def beamDynToHawc2(BD_mainfile, BD_bladefile, H2_htcfile=None, H2_stfile=None, bodyname=None, A=None, E=None, G=None, theta_p_in=None, FPM=False, verbose=False):
    """ 
    
     FPM: fully populated matrix, if True, use the FPM format of hawc2
    """
    # --- Read BeamDyn files
    if isinstance(BD_mainfile, str):
        BD_mainfile = FASTInputFile(BD_mainfile)
    if isinstance(BD_bladefile, str):
        BD_bladefile = FASTInputFile(BD_bladefile)
    bdLine = BD_mainfile.toDataFrame()
    bd     = BD_bladefile.toDataFrame()

    # --- Extract relevant info
    prop  = bd['BeamProperties']
    kp_x  = bdLine['kp_xr_[m]'].values
    kp_y  = bdLine['kp_yr_[m]'].values
    kp_z  = bdLine['kp_zr_[m]'].values
    twist = bdLine['initial_twist_[deg]'].values*np.pi/180 # BeamDyn convention
    r_bar = prop['Span'].values

    K = np.zeros((6,6),dtype='object')
    M = np.zeros((6,6),dtype='object')
    for i in np.arange(6):
        for j in np.arange(6):
            K[i,j]=prop['K{}{}'.format(i+1,j+1)].values
            M[i,j]=prop['M{}{}'.format(i+1,j+1)].values

    # Map 6x6 data to "beam" data
    # NOTE: theta_* are in [rad]
    EA, EIx, EIy, kxsGA, kysGA, GKt, x_C, y_C, x_S, y_S, theta_p, theta_s = K66toProps(K, theta_p_in)
    m, Ixi, Iyi, Ip, x_G, y_G, theta_i = M66toProps(M)
#     print('kxGA    {:e}'.format(np.mean(kxsGA)))
#     print('kyGA    {:e}'.format(np.mean(kysGA)))
#     print('EA      {:e}'.format(np.mean(EA)))
#     print('EIx     {:e}'.format(np.mean(EIx)))
#     print('EIy     {:e}'.format(np.mean(EIy)))
#     print('GKt     {:e}'.format(np.mean(GKt)))
#     print('xC    ',np.mean(x_C))
#     print('yC    ',np.mean(y_C))
#     print('xS    ',np.mean(x_S))
#     print('yS    ',np.mean(y_S))
#     print('thetap',np.mean(theta_p))
#     print('thetas',np.mean(theta_s))
#     print('m     ',np.mean(m))
#     print('Ixi   ',np.mean(Ixi))
#     print('Iyi   ',np.mean(Iyi))
#     print('Ip    ',np.mean(Ip))
#     print('x_G   ',np.mean(x_G))
#     print('y_G   ',np.mean(y_G))
#     print('thetai',np.mean(theta_i))

    # Convert to Hawc2 system
    if FPM:
        dfMeanLine , dfStructure = beamDyn2Hawc2FPM_raw(r_bar,
                kp_x, kp_y, kp_z, twist,  # BeamDyn convention, twist around -z [in rad]
                m, Ixi, Iyi, x_G, y_G, theta_i,  # theta_i/p around z (in rad)
                x_C, y_C, theta_p, K)

    else:

        dfMeanLine , dfStructure = beamDyn2Hawc2_raw(r_bar,
                kp_x, kp_y, kp_z, twist, 
                m, Ixi, Iyi, x_G, y_G, theta_i,
                EA, EIx, EIy, GKt, kxsGA, kysGA, x_C, y_C, theta_p, x_S, y_S, theta_s, 
                A=A, E=E, G=G)

    # --- Rewrite st file
    if H2_stfile is not None:
        try:
            os.makedirs(os.path.dirname(H2_stfile))
        except:
            pass
        if verbose: 
            print('Writing:   ',H2_stfile)
        with open(H2_stfile, 'w') as f:
            f.write('%i ; number of sets, Nset\n' % 1)
            f.write('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            f.write('#%i ; set number\n' % 1)
            if FPM:
                cols=['r','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
            else:
                cols=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]', 'x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
            f.write('\t'.join(['{:20s}'.format(s) for s in cols])+'\n')
            f.write('$%i %i\n' % (1, dfStructure.shape[0]))
            f.write('\n'.join('\t'.join('%19.13e' %x for x in y) for y in dfStructure.values))

    # --- Rewrite htc file
    if H2_htcfile is not None:
        def readToMarker(lines, marker, i, nMax=None, noException=False):
            l_sel=[]
            if nMax is None: nMax=len(lines)
            while i<nMax:
                line=lines[i]
                if line.replace(' ','').lower().find(marker)>=0:
                    break
                l_sel.append(line.strip())
                i+=1
            if line.strip().replace(' ','').lower().find(marker)<0:
                if noException:
                    return None, None, None
                else:
                    raise Exception('Marker not found '+ marker)
            return l_sel, line, i

        with open(H2_htcfile, 'r') as f:
            lines_in = f.readlines()
        lines_out = []
        bodyNotFound=True
        iBodyEnd=0
        nBodies=0
        while bodyNotFound and nBodies<10:
            _, line, iBodyStart = readToMarker(lines_in, 'beginmain_body',iBodyEnd)
            _, line, iBodyEnd = readToMarker(lines_in, 'endmain_body', iBodyStart)
            _, line, iBody = readToMarker(lines_in, 'name'+bodyname, iBodyStart, iBodyEnd, True)
            nBodies+=1
            if line is None:
                iBody=-1
            else:
                #print('Body {} found between lines {} and {} '.format(bodyname, iBodyStart+1, iBodyEnd+1))
                bodyNotFound=False
        if nBodies>=10:
            raise Exception('Body {} not found in file'.format(bodyname))

        _, line, iC2Start = readToMarker(lines_in, 'beginc2_def', iBodyStart, iBodyEnd)
        _, line, iC2End   = readToMarker(lines_in, 'endc2_def'  , iC2Start, iBodyEnd)

        _, line, iTIStart = readToMarker(lines_in, 'begintimoschenko_input', iBodyStart, iBodyEnd)
        _, line, iTIEnd   = readToMarker(lines_in, 'endtimoschenko_input'  , iTIStart, iBodyEnd)


        simdir        = os.path.dirname(H2_htcfile)
        H2_stfile_rel = os.path.relpath(H2_stfile, simdir)

        lines_out  = lines_in[:iTIStart+1]
        lines_out += ['      filename {};\n'.format(H2_stfile_rel)]
        if FPM:
            lines_out += ['      FPM 1 ;\n']
        #    lines_out += ['      FPM 0 ;\n']
        lines_out += ['      set 1 1 ;\n']
        lines_out += lines_in[iTIEnd:iC2Start+1]
        lines_out += ['      nsec {} ;\n'.format(dfMeanLine.shape[0])]
        for i, row in dfMeanLine.iterrows():
            lines_out += ['      sec {:4d}\t{:13.6e}\t{:13.6e}\t{:13.6e}\t{:13.6e};\n'.format(i+1, row['x_[m]'],row['y_[m]'],row['z_[m]'],row['twist_[deg]'])]
        lines_out += lines_in[iC2End:]


        if verbose: 
            print('ReWriting: ',H2_htcfile)
        with open(H2_htcfile, 'w') as f:
            f.write(''.join(lines_out))

    return dfMeanLine, dfStructure


def beamDyn2Hawc2FPM_raw(r_bar, kp_x, kp_y, kp_z, twist,
        m, Ixi, Iyi, x_G, y_G, theta_i, 
        x_C, y_C, theta_p,  
        K):
    """
    NOTE: all angles are in radians
    
    """
    import scipy.linalg
    # --- BeamDyn to Hawc2 Structural data
#     import pdb; pdb.set_trace()
    # Hawc2 = BeamDyn
    x_cg    = -y_G
    y_cg    = x_G
    x_e     = -y_C
    y_e     = x_C
    pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    if np.all(np.abs(m)<1e-16):
        ri_y    = m*0
        ri_x    = m*0
    else:
        ri_y    = np.sqrt(Ixi/m)    # [m]
        ri_x    = np.sqrt(Iyi/m)    # [m]
    # Curvilinear position of keypoints (only used to get max radius...)
    dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    r_p= np.concatenate(([0],np.cumsum(dr)))
    r=r_bar * r_p[-1]

    RotMat=np.array([  # From Hawc2 to BeamDyn
            [0 ,1,0],
            [-1,0,0],
            [0,0,1]])
    RR= scipy.linalg.block_diag(RotMat,RotMat)

    nSpan = len(K[0,0])
    KH2=np.zeros((6,6,nSpan))
    for iSpan in np.arange(nSpan):
        Kbd = np.zeros((6,6))
        for i in np.arange(6):
            for j in np.arange(6):
                Kbd[i,j] = K[i,j][iSpan]
        Kh2 = (RR.T).dot(Kbd).dot(RR)
        for i in np.arange(6):
            for j in np.arange(6):
                KH2[i,j][iSpan]=Kh2[i,j]

    K11 = KH2[0,0]
    K22 = KH2[1,1]
    K33 = KH2[2,2]
    K44 = KH2[3,3]
    K55 = KH2[4,4]
    K66 = KH2[5,5]

    K12 = KH2[0,1]
    K13 = KH2[0,2]
    K14 = KH2[0,3]
    K15 = KH2[0,4]
    K16 = KH2[0,5]
    K23 = KH2[1,2]
    K24 = KH2[1,3]
    K25 = KH2[1,4]
    K26 = KH2[1,5]
    K34 = KH2[2,3]
    K35 = KH2[2,4]
    K36 = KH2[2,5]
    K44 = KH2[3,3]
    K45 = KH2[3,4]
    K46 = KH2[3,5]
    K55 = KH2[4,4]
    K56 = KH2[4,5]


    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
    data = np.column_stack((r, m, x_cg, y_cg, ri_x, ri_y, pitch, x_e, y_e, K11,K12,K13,K14,K15,K16,K22,K23,K24,K25,K26,K33,K34,K35,K36,K44,K45,K46,K55,K56,K66))
    dfStructure = pd.DataFrame(data=data, columns=columns)

#     #      Siemens z ->  BeamDyn z
#     #      Siemens x -> -BeamDyn y
#     #      Siemens y ->  BeamDyn x
#     KSiemens = np.zeros((nSpan,6,6))
#     for i in np.arange(6):
#         for j in np.arange(6):
#             if j>=i:
#                 key='d{}{}'.format(i+1,j+1)
#             else:
#                 key='d{}{}'.format(j+1,i+1)
#             KSiemens[:,i,j] = sd[key]
# 
#     for i in np.arange(len(bp['span'])):
#         K = bp['K'][i]*0
#         Ks= KSiemens[i]
#         K = RR.dot(Ks).dot(RR.T)
#         bp['K'][i] = K



    # --- BeamDyn to Hawc2 Reference axis
    X_H2     = -kp_y
    Y_H2     = kp_x
    Z_H2     = kp_z
    twist_H2 = - twist*180/np.pi # - [deg]
    columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure



def beamDyn2Hawc2_raw(r_bar, kp_x, kp_y, kp_z, twist,
        m, Ixi, Iyi, x_G, y_G, theta_i, 
        EA, EIx, EIy, GKt, kxsGA, kysGA, x_C, y_C, theta_p, x_S, y_S, theta_s, 
        A=None, E=None, G=None):
    """ 
    NOTE: all angles are in radians
    """
    # --- BeamDyn to Hawc2 Structural data
    if A is None: A = np.ones(x_G.shape)
    if E is None: E = EA/A
    if G is None: G = E/2/(1+0.3) # Young modulus
#     import pdb; pdb.set_trace()
    # Hawc2 = BeamDyn
    x_cg    = -y_G
    y_cg    = x_G
    x_sh    = -y_S
    y_sh    = x_S
    x_e     = -y_C
    y_e     = x_C
    I_y     = EIx/E            # [m^4] Hawc2 Iy is wrt to principal bending ye axis
    I_x     = EIy/E            # [m^4] Hawc2 Ix is wrt to principal bending xe axis
    I_p     = GKt/G            # [m^4]
    k_y     = kxsGA/(G*A)
    k_x     = kysGA/(G*A)
    pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    if np.all(np.abs(m)<1e-16):
        ri_y    = m*0
        ri_x    = m*0
    else:
        ri_y    = np.sqrt(Ixi/m)    # [m]
        ri_x    = np.sqrt(Iyi/m)    # [m]
    # Curvilinear position of keypoints (only used to get max radius...)
    dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    r_p= np.concatenate(([0],np.cumsum(dr)))
    r=r_bar * r_p[-1]

    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
    data = np.column_stack((r,m,x_cg,y_cg,ri_x, ri_y, x_sh, y_sh, E, G, I_x, I_y, I_p, k_x, k_y, A, pitch, x_e, y_e))
    dfStructure = pd.DataFrame(data=data, columns=columns)

    # --- BeamDyn to Hawc2 Reference axis
    X_H2     = -kp_y
    Y_H2     = kp_x
    Z_H2     = kp_z
    twist_H2 = - twist*180/np.pi # -[deg]
    columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure




# --------------------------------------------------------------------------------}
# --- Functions for 6x6 matrices 
# --------------------------------------------------------------------------------{
def M66toProps(M, convention='BeamDyn'):
    """ 
    Convert mass properties of a 6x6 section to beam properties
    This assumes that the axial and bending loads are decoupled.

    INPUTS:
     - M : 6x6 array of mass elements. Each element may be an array (e.g. for all spanwise values)
     OUTPUTS:
      - m: section mass
      - Ixx, Iyy, Ixy: area moment of inertia
      - x_G, y_G
    """
    M11=M[0,0]
    M44=M[3,3]
    M55=M[4,4]
    M66=M[5,5]
    M16=M[0,5]
    M26=M[1,5]
    M45=M[3,4]

    m=M11
    if convention=='BeamDyn':
        if np.all(np.abs(m)<1e-16):
            Ixi = Iyi = Ipp = x_G = y_G = theta_i = m*0
            return m, Ixi, Iyi, Ipp, x_G, y_G, theta_i
        y_G= -M16/m
        x_G=  M26/m
        # sanity
        np.testing.assert_array_almost_equal([M[0,3],M[0,4]],[0*M[0,3],0*M[0,3]])
        np.testing.assert_array_almost_equal([M[1,3],M[1,4]],[0*M[0,3],0*M[0,3]])
        np.testing.assert_array_almost_equal([M[3,5],M[4,5]],[0*M[0,3],0*M[0,3]])
        
        Ixx =  M44-m*y_G**2
        Iyy =  M55-m*x_G**2
        Ixy = -M45-m*x_G*y_G
        Ipp =  M66 -m*(x_G**2 + y_G**2)

        if np.all(np.abs(Ixy)<1e-16):
            # NOTE: Assumes theta_i ==0
            #print('>>> Assume theta_i 0')
            Ixi     = Ixx
            Iyi     = Iyy
            theta_i = Ixx*0
        else:
            #print('>>> Minimize theta_i')
            Ixi = np.zeros(Ixx.shape)
            Iyi = np.zeros(Ixx.shape)
            theta_i = np.zeros(Ixx.shape)
            for i, (hxx, hyy, hxy) in enumerate(zip(Ixx,Iyy,Ixy)):
                Ixi[i],Iyi[i],theta_i[i] = solvexytheta(hxx,hyy,hxy)

                MM2= MM(m[i],Ixi[i],Iyi[i],Ipp[i],x_G[i],y_G[i],theta_i[i])

                np.testing.assert_allclose(MM2[3,3], M[3,3][i], rtol=1e-3)
                np.testing.assert_allclose(MM2[4,4], M[4,4][i], rtol=1e-3)
                np.testing.assert_allclose(MM2[5,5], M[5,5][i], rtol=1e-3)
                np.testing.assert_allclose(MM2[3,4], M[3,4][i], rtol=1e-3)

        np.testing.assert_array_almost_equal(Ipp, Ixx+Iyy, 2)
        np.testing.assert_array_almost_equal(Ipp, Ixi+Iyi, 2)
         
    else:
        raise NotImplementedError()

    return m, Ixi, Iyi, Ipp, x_G, y_G, theta_i

def solvexytheta(Hxx,Hyy,Hxy):
    """ 
     Solve for a system of three unknown, used to get:
       - EIy, EIx and thetap given Hxx,Hyy,Hxy
       - kxs*GA, kys*GA and thetas given Kxx,Kyy,Kxy
       - I_x, I_y and theta_is given Ixx,Iyy,Ixy
    """
    from scipy.optimize import fsolve
    def residual(x):
        EI_x, EI_y, theta_p =x
        res=np.array([
            Hxx - EI_x*np.cos(theta_p)**2 - EI_y*np.sin(theta_p)**2 ,
            Hyy - EI_x*np.sin(theta_p)**2 - EI_y*np.cos(theta_p)**2,
            Hxy - (EI_y-EI_x)*np.sin(theta_p)*np.cos(theta_p)]
                ).astype(float)
        return res
    x0 = [Hxx,Hyy,0]
    x_opt   = fsolve(residual, x0)
    EI_x,EI_y,theta_p = x_opt
    theta_p = np.mod(theta_p,2*np.pi)
    return EI_x, EI_y, theta_p



def K66toProps(K, theta_p_in=None, convention='BeamDyn'):
    """ 
    Convert stiffness properties of a 6x6 section to beam properties
    This assumes that the axial and bending loads are decoupled.

    INPUTS:
     - K : 6x6 array of stiffness elements. Each element may be an array (e.g. for all spanwise values)
    INPUTS OPTIONAL:
     - theta_p_in : angle from section to principal axis [rad], positive around z
     - convention : to change coordinate systems in the future
    OUTPUTS:
     - EA, EIx, EIy: axial and bending stiffnesses
     - kxGA, kyGA, GKt: shear and torsional stiffness
     - xC,yC : centroid
     - xS,yS : shear center
     - theta_p, theta_s: angle to principal axes and shear axes [rad]
    """
    K11=K[0,0]
    K22=K[1,1]
    K33=K[2,2]
    K44=K[3,3]
    K55=K[4,4]
    K66=K[5,5]

    K12=K[0,1]
    K16=K[0,5]
    K26=K[1,5]
    K34=K[2,3]
    K35=K[2,4]
    K45=K[3,4]

    if convention=='BeamDyn':
        # --- EA, EI, centroid, principal axes
        EA =  K33
        yC =  K34/EA
        xC = -K35/EA
        Hxx=  K44-EA*yC**2
        Hyy=  K55-EA*xC**2 # NOTE: xC fixed
        Hxy= -K45-EA*xC*yC # NOTE: sign changed

        if theta_p_in is not None:
            theta_p=theta_p_in
            print('>>> theta_p given')
            C2=np.cos(theta_p)**2
            S2=np.sin(theta_p)**2
            C4=np.cos(theta_p)**4
            S4=np.sin(theta_p)**4
            EIxp = (Hxx*C2 - Hyy*S2)/(C4-S4)
            EIyp = (Hxx*S2 - Hyy*C2)/(S4-C4)
            Hxyb = (EIyp-EIxp)*np.sin(theta_p)*np.cos(theta_p)

            bNZ=np.logical_and(Hxy!=0, Hxyb!=0)
            np.testing.assert_allclose(Hxy[bNZ], Hxyb[bNZ], rtol=1e-3)
            np.testing.assert_allclose(EIxp+EIyp, Hxx+Hyy, rtol=1e-3)

        else:
            if np.all(np.abs(Hxy)<1e-16):
                #print('>>>> assume theta_p=0')
                # NOTE: Assumes theta_p ==0
                EIxp = Hxx
                EIyp = Hyy
                theta_p=0*EA
            else:
                #print('>>> Minimization for theta_p')
                EIxp= np.zeros(Hxx.shape)
                EIyp= np.zeros(Hxx.shape)
                theta_p = np.zeros(Hxx.shape)
                for i, (hxx, hyy, hxy) in enumerate(zip(Hxx,Hyy,Hxy)):
                    EIxp[i],EIyp[i],theta_p[i] = solvexytheta(hxx,hyy,hxy)

            theta_p[theta_p>np.pi]=theta_p[theta_p>np.pi]-2*np.pi

        # --- Torsion, shear terms, shear center
        Kxx =  K11
        Kxy = -K12
        Kyy =  K22
        yS  = (Kyy*K16+Kxy*K26)/(-Kyy*Kxx + Kxy**2)
        xS  = (Kxy*K16+Kxx*K26)/( Kyy*Kxx - Kxy**2)
        GKt = K66 - Kxx*yS**2 -2*Kxy*xS*yS - Kyy*xS**2
        if np.all(np.abs(Kxy)<1e-16):
            # Assumes theta_s=0
            kxsGA = Kxx # Kxx = kxs*GA
            kysGA = Kyy
            theta_s=0*EA
        else:
            kxsGA = np.zeros(Kxx.shape)
            kysGA = np.zeros(Kxx.shape)
            theta_s = np.zeros(Hxx.shape)
            for i, (kxx, kyy, kxy) in enumerate(zip(Kxx,Kyy,Kxy)):
                kxsGA[i],kysGA[i],theta_s[i] = solvexytheta(kxx,kyy,kxy)

        theta_s[theta_s>np.pi]=theta_s[theta_s>np.pi]-2*np.pi


        # sanity checks
        KK2= KK(EA, EIxp, EIyp, GKt, EA*0+1, kxsGA, kysGA, xC, yC, theta_p, xS, yS, theta_s)
        np.testing.assert_allclose(KK2[0,0], K[0,0], rtol=1e-2)
        np.testing.assert_allclose(KK2[1,1], K[1,1], rtol=1e-2)
        np.testing.assert_allclose(KK2[2,2], K[2,2], rtol=1e-2)
        np.testing.assert_allclose(KK2[3,3], K[3,3], rtol=1e-1)
#         np.testing.assert_allclose(KK2[4,4], K[4,4], rtol=1e-2)
        np.testing.assert_allclose(KK2[5,5], K[5,5], rtol=1e-1)
        np.testing.assert_allclose(KK2[2,3], K[2,3], rtol=1e-2)
        np.testing.assert_allclose(KK2[2,4], K[2,4], rtol=1e-2)

        np.testing.assert_allclose(K16, -Kxx*yS-Kxy*xS)
#         np.testing.assert_allclose(KK2[0,5], K[0,5],rtol=1e-3)
#         np.testing.assert_allclose(KK2[1,5], K[1,5],rtol=5e-2) # Kxy harder to get
        #np.testing.assert_allclose(KK2[3,4], K[3,4]) # <<< hard to match

    else:
        raise NotImplementedError()

    return EA, EIxp, EIyp, kxsGA, kysGA, GKt, xC, yC, xS, yS, theta_p, theta_s


if __name__=='__main__':
    np.set_printoptions(linewidth=300)

    # --- BeamDyn 2 Hawc 2
    BD_mainfile  = 'solid_beam_BeamDyn.dat'
    BD_bladefile = '../solid_beam_BeamDyn_Blade.dat'
    H2_htcfile_old  = './_template.htc'
    H2_htcfile_new  = './solid_beam_hawc2.htc'
    H2_stfile    = './solid_beam_st.dat'

    from shutil import copyfile
    copyfile(H2_htcfile_old, H2_htcfile_new)

    beamDyn2Hawc2(BD_mainfile, BD_bladefile, H2_htcfile_new, H2_stfile, 'beam_1', A=None, E=None, G=None)


