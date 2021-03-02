#############################################################
#	Extract OpenFAST Matrices from linearization files  #
#	Authors: Srinivasa B. Ramisetti    	            #
#	Created:   01-July-2020		  	            #
#	E-mail: ramisettisrinivas@yahoo.com		    #
#	Web:	http://ramisetti.github.io		    #
#############################################################
#!/usr/bin/env python

import os, distutils
import collections
from itertools import islice
import numpy as np
import re

def _isbool(str):
    flag=0
    if str=='T' or str=='True':
        flag=1
    return flag


def findBladeTriplets_EDstate(rotFrame,Desc):
    # descriptions of ED states contain the blade number twice. This is the old
    # method of comparing strings.
    #   frame:
    NTriplets = 0;                  # first initialize to zero 
    Triplets = [];
    for i in range(len(rotFrame)):  # loop through all active (enabled) DOFs 
        if rotFrame[i]:  # this is a state in the rotating frame 

            col = int(Desc[i].find('blade'))              # find the starting index of the string 'blade'
            if col!=-1:                             # true if the Desc{i} contains the string 'blade' 
                k = int(Desc[i][col+6]);          # save the blade number for the initial blade
                Tmp = -1*np.ones(3);
                #Tmp = np.zeros(3);                        # first initialize to zero 
                Tmp[k-1] = int(i);                              # save the index for the initial blade 

                # find the other match values

                for j in range((i+1),len(rotFrame)):           # loop through all remaining active (enabled) DOFs 
                    if Desc[j][:col]==Desc[i][:col]:      # true if we have the same state from a different blade 
                        k = int(Desc[j][col+6]);  # save the blade numbers for the remaining blades 

                        Tmp[k-1] = int(j);                      # save the indices for the remaining blades
                        TmpTmp=Tmp+1
                        if ( (Tmp>-1).all() ):
                        #if ( np.all(TmpTmp+1) ):                  # true if all the elements of Tmp are nonzero; thus, we found a triplet of rotating indices 
                            NTriplets = NTriplets + 1;   # this  is  the number  of  state triplets in the rotating frame 
                            Triplets.append(Tmp); # these are the indices for state triplets in the rotating frame
                            break
    
    return Triplets,NTriplets;


def findBladeTriplets(rotFrame,Desc):

    # Find the number of, and indices for, triplets in the rotating frame:
    chkStr = ['[Bb]lade \d', '[Bb]lade [Rr]oot \d', 'BD_\d', '[Bb]\d', '[Bb]lade\d', 'PitchBearing\d', '\d']

    NTriplets = 0;              # first initialize to zero
    Triplets = [];
    for i in range(len(rotFrame)):  # loop through inputs/outputs
    #for i in range(7):
        #if(i>=67):
            #print(rotFrame[i])
        if rotFrame[i]:          # this is in the rotating frame
            Tmp = -1*np.ones(3);
            foundTriplet = False;
            foundIt = False;
            for chk in chkStr:
                BldNoCol = re.search(chk,Desc[i]);
                if BldNoCol!=None:
                    foundIt = True;

                    Bldstr=BldNoCol.group()
                    # create another regular expression to find the
                    # exact match on a different blade:
                    strng = re.split(Bldstr,Desc[i],1); #this should return the strings before and after the match
                    #print(strng[1])

                    
                    FirstStr = strng[0] + Bldstr[:len(Bldstr)-1] + '.'
                    checkThisStr = FirstStr + strng[1]
                    
                    #we need to get rid of the special characters that
                    #may exist in Desc{}:
                    checkThisStr=checkThisStr.replace(')','\)').replace('(', '\(').replace('^','\^')

                    k = int(Bldstr[len(Bldstr)-1])
                    #print(Bldstr[len(Bldstr)-1], checkThisStr)
                    Tmp[k-1] = int(i);
                    #print(Tmp,k,i,len(rotFrame),foundIt)
                    #exit()
                    break

                #print(Tmp,j)
            
            # find the other match values
            if foundIt:
                for j in range((i+1),len(rotFrame)):           # loop through all remaining control inputs
                    #print(i, j)
                    if rotFrame[j]:                       # this is in the rotating frame
                        BldNoCol = re.search(checkThisStr, Desc[j])
                        if BldNoCol!=None:
                            Num = re.search(FirstStr,Desc[j]).group();
                            #Num = regexp(Desc{j},FirstStr,'match'); # match all but the blade number
                            k = int(Num[len(Num)-1]);
                            Tmp[k-1] = int(j);                             # save the indices for the remaining blades
                            #TmpTmp=Tmp+1
                            #print(BldNoCol.group(),i,j,k)
                            if ( (Tmp>-1).all() ):                         # true if all the elements of Tmp are nonzero; thus, we found a triplet of rotating indices
                                foundTriplet = True;
                                
                                Triplets.append(Tmp);        # these are the indices for control input triplets in the rotating frame

                                NTriplets = NTriplets + 1;          # this  is  the number  of  control input triplets in the rotating frame

                                # we'll set rotFrame to false so that we don't have to check the found channels again; also allows us to throw error if we have a rotating channel that doesn't have a unique match
                                for idx in Tmp:
                                    id=int(idx)
                                    rotFrame[id] = 0;
                                
                                break;
                
                if foundTriplet==False:
                    print('Rotating channel "', i, Desc[i], '" does not form a unique blade triplet. Blade(s) not found: ', np.array(np.where(Tmp == -1))+1 )
            else:
                print( 'Could not find blade number in rotating channel "', Desc[i], '".')
    #print(NTriplets)
    #print(Triplets)

    return Triplets, NTriplets

def reOrderByIdx_1D(arry,Indi):
    tmp=[None]*len(arry)
    j=0
    for i in Indi:
        tmp[i]=arry[j]
        j=j+1
    return tmp

def reOrderByIdx_2D(arry,Indi,Indj):
    #tmp=[[None]*arry.shape[0]]*arry.shape[1]
    tmp=np.empty((arry.shape[0], arry.shape[1])) 
    #print(arry.shape, len(Indi), len(Indj))
    if arry.shape[0]!= len(Indi):
        Indi=np.arange(0,arry.shape[0])
    if arry.shape[1]!= len(Indj):
        Indj=np.arange(0,arry.shape[1])

    p=0
    q=0
    for i in Indi:
        q=0
        for j in Indj:
            tmp[i][j]=arry[p][q]
            q=q+1
        p=p+1
    return tmp

def getStateOrderingIndx(matData):

    #StateOrderingIndx={}
    StateOrderingIndx = np.arange(0,matData['NumStates'])
    lastModName = '';
    checkEDstates = True;
    lastModOrd  = 0;
    mod_nDOFs   = 0;    # number of DOFs in each module
    sum_nDOFs2  = 0;    # running total of second-order DOFs
    sum_nDOFs1  = 0;    # running total of first-order DOFs
    indx_start  = 0;    # starting index of the modules
    
    for i in range(0,matData['NumStates']):

        tmp=(matData['DescStates'][i]); # name of the module whose states we are looking at
        modName = tmp.split(' ')[0]

        #print('indx ', i)
        # ED has the blade number in the state description twice, so we
        # have to check the strings differently. We'll check here if a  
        # different module is used for the blade DOFs:
        if modName[:2] =='BD' or modName[:3]=='MBD':
            checkEDstates = False;

        if lastModName!=modName:
            # this is the start of a new set of DOFs, so we'll set the
            # indices for the last matrix
            if lastModOrd == 2:
                mod_nDOFs = int(mod_nDOFs/2);
                #print('mmmm ', mod_nDOFs, indx_start, matData['ndof2'], i)
                StateOrderingIndx[  indx_start           :(indx_start+mod_nDOFs)] = sum_nDOFs2 +                 np.arange(0,mod_nDOFs); # q2 starts at 1
                StateOrderingIndx[ (indx_start+mod_nDOFs):(i)]                  = sum_nDOFs2 + matData['ndof2'] + np.arange(0,mod_nDOFs); # q2_dot starts at matData.ndof2 + 1

                sum_nDOFs2 = sum_nDOFs2 + mod_nDOFs;
            else:
                if indx_start < mod_nDOFs:
                    StateOrderingIndx[indx_start:(indx_start+mod_nDOFs)] = sum_nDOFs1 + matData['NumStates2'] + np.arange(0,mod_nDOFs); # q1 starts at matData.NumStates2 + 1
                
                sum_nDOFs1 = sum_nDOFs1 + mod_nDOFs;
            
            # reset for a new module
            mod_nDOFs = 0;
            
            indx_start = i; #start of this module
            lastModName = modName;
            lastModOrd  = matData['StateDerivOrder'][i];
            
        mod_nDOFs = mod_nDOFs+1;
        

    # repeat for the last module found:
    if lastModOrd == 2:
        mod_nDOFs = int(mod_nDOFs/2);
        #print(mod_nDOFs,indx_start)
        StateOrderingIndx[indx_start:(indx_start+mod_nDOFs)] = sum_nDOFs2 + np.arange(0,mod_nDOFs); # q2 starts at 1
        StateOrderingIndx[(indx_start+mod_nDOFs):matData['NumStates']]        = sum_nDOFs2 + matData['ndof2'] + np.arange(0,mod_nDOFs); # q2_dot starts at matData.ndof2 + 1
    else:
        StateOrderingIndx[indx_start:(indx_start+mod_nDOFs)] = sum_nDOFs1 + matData['NumStates2'] + np.arange(0,mod_nDOFs); # q1 starts at matData.NumStates2 + 1

    #print(StateOrderingIndx)
    return StateOrderingIndx, checkEDstates


def readLinTable(f,n):
    tmp={}
    op=[]
    RF=[]
    DO=[]
    desc=[]
    for line in islice(f,n):
        op.append(np.float(line.split()[1]))
        RF.append(_isbool(line.split()[2]))
        DO.append(np.int(line.split()[3]))
        desc.append(" ".join(line.split()[4:]))
        
    tmp['x_op']=op
    tmp['x_rotFrame']=RF
    tmp['x_DerivOrder']=DO
    tmp['x_desc']=desc
    return tmp

def readFASTMatrix(f):
    name=""
    m=0
    tmp=[]
    for line in islice(f,1):
        # copy matrix name to tmp_name
        tmp_name=line.strip().split(':')[0]
        # get matrix dimensions if matrix exists
        if tmp_name != "":
            m=np.int(line.strip().split(':')[1].split('x')[0])
            n=np.int(line.strip().split(':')[1].split('x')[1])

    # copy matrix into tmp list
    if m!=0:
        name=tmp_name
        for line in islice(f,m):
            tmp.append([np.float(num) for num in line.split()])
        tmp=np.array(tmp)
    return name,tmp
    
def ReadFASTLinear(filename):
    with open(filename) as f:
        info = {}
        data = {}
        SetOfMatrices = 1
        info['name'] = os.path.splitext(os.path.basename(filename))[0]
        try:
            header = [f.readline() for _ in range(17)]
            info['fast_version'] = header[1].strip()
            info['modules'] = header[2].strip()
            info['description'] = header[4].strip()

            ln=7;
            #data = np.array([line.split() for line in f.readlines()]).astype(np.float)
            data['t']=np.float(header[ln].split(':')[1].strip().split(' ')[0])
            data['RotSpeed']=np.float(header[ln+1].split(':')[1].strip().split(' ')[0])
            data['Azimuth']=np.float(header[ln+2].split(':')[1].strip().split(' ')[0])
            data['WindSpeed']=np.float(header[ln+3].split(':')[1].strip().split(' ')[0])
            data['n_x']=np.float(header[ln+4].split(':')[1].strip().split(' ')[0])
            data['n_xd']=np.float(header[ln+5].split(':')[1].strip().split(' ')[0])
            data['n_z']=np.float(header[ln+6].split(':')[1].strip().split(' ')[0])
            data['n_u']=np.float(header[ln+7].split(':')[1].strip().split(' ')[0])
            data['n_y']=np.float(header[ln+8].split(':')[1].strip().split(' ')[0])
            
            data['Azimuth']=np.mod(data['Azimuth'],2.0*np.pi)

            if header[ln+9].split('?')[1].strip()=='Yes':
                SetOfMatrices=2

            # skip next three lines
            for line in islice(f,4):
                pass
            
            #ln=ln+9+3

            #header = [f.readline() for _ in range(ln,3)]
            
            if data['n_x'] > 0:
                temp=readLinTable(f,int(data['n_x']))
                data['x_op']=temp['x_op']
                data['x_rotFrame']=temp['x_rotFrame']
                data['x_DerivOrder']=temp['x_DerivOrder']
                data['x_desc']=temp['x_desc']

                # skip next three lines
                for line in islice(f,4):
                    pass

                temp=readLinTable(f,int(data['n_x']))
                data['xdot_op']=temp['x_op']
                data['xdot_desc']=temp['x_desc']

                if data['x_DerivOrder'][1] == 0: # this is an older file without derivOrder columns
                    # (these are second-order states)
                    data['x_DerivOrder']=2.0
                    data['n_x2'] = data['n_x'] 
                else:
                    #(number of second-order states)
                    data['n_x2'] = sum(1 for i in data['x_DerivOrder'] if i == 2)
            else:
                data['n_x2'] = 0;
            
            
            if data['n_xd'] > 0:
                # skip next three lines
                for line in islice(f,4):
                    pass
                temp=readLinTable(f,int(data['n_xd']))
                data['xd_op']=temp['x_op']
                data['xd_desc']=temp['x_desc']
            if data['n_z'] > 0:
                # skip next three lines
                for line in islice(f,4):
                    pass
                temp=readLinTable(f,int(data['n_z']))
                data['z_op']=temp['x_op']
                data['z_desc']=temp['x_desc']
            if data['n_u'] > 0:
                # skip next three lines
                for line in islice(f,4):
                    pass
                temp=readLinTable(f,int(data['n_u']))
                data['u_op']=temp['x_op']
                data['u_desc']=temp['x_desc']
                data['u_rotFrame']=temp['x_rotFrame']
            if data['n_y'] > 0:
                # skip next three lines
                for line in islice(f,4):
                    pass
                temp=readLinTable(f,int(data['n_y']))
                data['y_op']=temp['x_op']
                data['y_desc']=temp['x_desc']
                data['y_rotFrame']=temp['x_rotFrame']

            # skip next one line
            for line in islice(f,4):
                pass

            mat_names=[]
            while True:
                name,mat=readFASTMatrix(f)
                if not name:
                    break;
                mat_names.append(name)
                data[name]=mat

            return data, info

        except (ValueError, AssertionError):
            raise


def get_Mats(FileNames):
    matData={}
    matData['NAzimStep']= len(FileNames);
    data=[None]*matData['NAzimStep']
    NAzimStep=matData['NAzimStep']
    #print('NAzimStep : ', matData['NAzimStep'])

    # Input data from linearization files
    data[NAzimStep-1],info     = ReadFASTLinear(FileNames[NAzimStep-1]);
    matData['NumStates']       = int(data[NAzimStep-1]['n_x']);
    matData['NumStates2']      = int(data[NAzimStep-1]['n_x2']);

    #return data, matData
    matData['ndof1']           = int(matData['NumStates'] - matData['NumStates2']); # number of first-order states = number of first-order DOFs
    matData['ndof2']           = int(data[NAzimStep-1]['n_x2'] / 2); # half the number of second-order states = number of second-order DOFs
    #% matData['ndof']            = matData['ndof2'] + matData['ndof1']; #half the number of second-order states plus the number of first-order states (i.e., states that aren't derivatives)

    matData['NumInputs']       = int(data[NAzimStep-1]['n_u']);
    matData['NumOutputs']      = int(data[NAzimStep-1]['n_y']);

    # allocate space for these variables
    
    matData['Azimuth']   = np.zeros(NAzimStep);
    matData['Omega']     = np.zeros(NAzimStep);
    matData['OmegaDot']  = np.zeros(NAzimStep);

    matData['WindSpeed'] = np.zeros(NAzimStep);

    if matData['NumStates'] > 0:
        matData['DescStates']      = data[NAzimStep-1]['x_desc'];
        matData['StateDerivOrder'] = data[NAzimStep-1]['x_DerivOrder'];
        matData['xdop']            = np.zeros((matData['NumStates'], NAzimStep))
        matData['xop']             = np.zeros((matData['NumStates'], NAzimStep))
        matData['A']               = np.zeros((matData['NumStates'], matData['NumStates'], NAzimStep))

    if matData['NumInputs'] > 0:
        matData['DescCntrlInpt'] = data[NAzimStep-1]['u_desc'];    
        matData['u_op']             = np.zeros((matData['NumInputs'],NAzimStep))
        if matData['NumStates']>0:
            matData['B'] = np.zeros((matData['NumStates'], matData['NumInputs'],NAzimStep))

    if matData['NumOutputs'] > 0:
        matData['DescOutput']    = data[NAzimStep-1]['y_desc']
        matData['y_op']             = np.zeros((matData['NumOutputs'],NAzimStep))

        if matData['NumStates'] > 0:
            matData['C']         = np.zeros((matData['NumOutputs'], matData['NumStates'], NAzimStep))
        if matData['NumInputs'] > 0:
            matData['D']         = np.zeros((matData['NumOutputs'], matData['NumInputs'], NAzimStep))

    # Reorder state matrices so that they follow the {q2, q2_dot, q1}
    # format that is assumed in the MBC3 equations.
    if matData['NumStates'] > 0:
        # keep StateOrderingIndx for applying inverse of MBC3 later
        # (to visualize mode shapes)
        matData['StateOrderingIndx'],checkEDstates    = getStateOrderingIndx(matData)

        sortedIndx=matData['StateOrderingIndx']
        #print(sortedIndx)
        x_rotFrame= reOrderByIdx_1D(data[NAzimStep-1]['x_rotFrame'],sortedIndx)
        matData['DescStates'] = reOrderByIdx_1D(data[NAzimStep-1]['x_desc'],sortedIndx)
        matData['StateDerivOrder'] = reOrderByIdx_1D(data[NAzimStep-1]['x_DerivOrder'],sortedIndx)

    for iFile in np.arange(0,NAzimStep):
        
        data[iFile],info     = ReadFASTLinear(FileNames[iFile]);
        matData['Omega'][iFile]   = data[iFile]['RotSpeed'];
        matData['Azimuth'][iFile] = data[iFile]['Azimuth']*180/np.pi;

        if 'WindSpeed' in matData:
            if 'WindSpeed' in data[iFile]:
                matData['WindSpeed'][iFile] = data[iFile]['WindSpeed']
            else:
                del matData['WindSpeed'];
    
        if 'A' in data[iFile]:
            matData['A'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['A'],sortedIndx,sortedIndx)
            #print('size of matData[A] for file ', iFile, ' is :', matData['A'].shape)
        if 'B' in data[iFile]:
            matData['B'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['B'],sortedIndx,np.arange(0,matData['NumStates']))
            #print('size of matData[B] for file ', iFile, ' is :', matData['B'].shape)
        if 'C' in data[iFile]:
            matData['C'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['C'],np.arange(0,matData['NumStates']),sortedIndx)
            #print('size of matData[C] for file ', iFile, ' is :', matData['C'].shape)
        if 'D' in data[iFile]:
            matData['D'][:,:,iFile]=reOrderByIdx_2D(data[iFile]['D'],sortedIndx,sortedIndx)
            #print('size of matData[D] for file ', iFile, ' is :', matData['D'].shape)

        if 'x_op' in data[iFile]:
            matData['xop'][:,iFile] = reOrderByIdx_1D(data[iFile]['x_op'],sortedIndx)
            #print('matData[xop] for file ', iFile,' is :',(matData['xop'][:,iFile]))
        if 'xdot_op' in data[iFile]:
            matData['xdop'][:,iFile] = reOrderByIdx_1D(data[iFile]['xdot_op'],sortedIndx)
            #print('matData[xdop] for file ', iFile,' is :',(matData['xdop'][:,iFile]))

        if 'u_op' in data[iFile]:
            matData['u_op'][:,iFile] = data[iFile]['u_op']

        if 'y_op' in data[iFile]:
            matData['y_op'][:,iFile] = data[iFile]['y_op']


    # Find the azimuth-averaged linearized 1st order state matrices:
    if 'A' in matData:
        matData['Avgxdop'] = np.mean(matData['xdop'],axis=1)
        matData['Avgxop']  = np.mean(matData['xop'], axis=1)
        matData['AvgA']    = np.mean(matData['A'],axis=2)

    #print(matData['AvgA'])
    if 'B' in matData:
        matData['AvgB']    = np.mean(matData['B'],axis=2)

    if 'C' in matData:
        matData['AvgC']    = np.mean(matData['C'],axis=2)

    if 'D' in matData:
        matData['AvgD']    = np.mean(matData['D'],axis=2)


    foundED  = True;
    for i in range(matData['ndof2']):
        # find the starting index of the string 'DOF_GeAz'
        if (matData['DescStates'][i].find('DOF_GeAz') != -1):
            matData['Omega']    = matData['xdop'][i,:]
            matData['OmegaDot'] = matData['xdop'][i+matData['ndof2'],:]
            foundED = True
            break

    for i in range(matData['ndof2']):
        # find the starting index of the string 'DOF_DrTr'
        if (matData['DescStates'][i].find('DOF_DrTr') != -1):
            matData['Omega']    = matData['Omega']    + matData['xdop'][i,:] #This always comes after DOF_GeAz so let's just add it here (it won't get written over later).
            matData['OmegaDot'] = matData['OmegaDot'] + matData['xdop'][i+matData['ndof2'],:]
            foundED = True
            break

    if not foundED:
        for i in range(matData['ndof2']):
            # find the starting index of the string 'Gearbox_Rot'
            if (matData['DescStates'][i].find('MBD Gearbox_Rot') != -1):
                matData['Omega']    = matData['xdop'][i,:]
                matData['OmegaDot'] = matData['xdop'][i+matData['ndof2'],:]
                break

    #print("\n".join(matData['DescStates']))
    #exit()
    # ----------- Find multi-blade coordinate (MBC) transformation indices ----

    # Find the indices for, state triplets in the rotating frame
    #  (note that we avoid the "first time derivative" states)
    if matData['ndof2'] > 0:
        if checkEDstates:
            matData['RotTripletIndicesStates2'], matData['n_RotTripletStates2'] = findBladeTriplets_EDstate(x_rotFrame[0:matData['ndof2']],matData['DescStates'][0:matData['ndof2']])
        else:
            matData['RotTripletIndicesStates2'], matData['n_RotTripletStates2'] = findBladeTriplets(x_rotFrame[0:matData['ndof2']],matData['DescStates'][0:matData['ndof2']])
    else:
        matData['RotTripletIndicesStates2'] = [];
        matData['n_RotTripletStates2'] = 0;

    if matData['ndof1'] > 0:
        matData['RotTripletIndicesStates1'], matData['n_RotTripletStates1'] = findBladeTriplets( x_rotFrame[matData['NumStates2']:] ,matData['DescStates'][matData['NumStates2']:] );
    else:
        matData['RotTripletIndicesStates1'] = [];
        matData['n_RotTripletStates1'] = 0;

    # Find the indices for control input triplets in the rotating frame:
    if matData['NumInputs'] > 0:
        matData['RotTripletIndicesCntrlInpt'], matData['n_RotTripletInputs'] = findBladeTriplets(data[0]['u_rotFrame'],matData['DescCntrlInpt'] );
    else:
        matData['RotTripletIndicesCntrlInpt'] = [];
        matData['n_RotTripletInputs'] = 0;

    # Find the indices for output measurement triplets in the rotating frame:
    if (matData['NumOutputs'] > 0 ):
        matData['RotTripletIndicesOutput'], matData['n_RotTripletOutputs'] = findBladeTriplets(data[0]['y_rotFrame'],matData['DescOutput'] );
    else:
        matData['RotTripletIndicesOutput'] = [];
        matData['n_RotTripletOutputs'] = 0;
        
    return matData, data

#matData,FAST_linData=get_Mats(FileNames)

#print("\n".join(matData['DescStates']))
                                      
# print(info)
# print(data['t'])
# print(data['RotSpeed'])
# print(data['Azimuth'])
# print(data['WindSpeed'])
# print(data['n_x'])
# print(data['n_xd'])
# print(data['n_z'])
# print(data['n_u'])
# print(data['n_y'])
# print(data['n_x2'])
# print(data['x_op'])
# print(data['x_desc'])
# print(data['x_rotFrame'])
# print(data['xdot_op'])
# print(data['xdot_desc'])
# print(data['u_op'])
# print(data['u_desc'])
# print(data['u_rotFrame'])
# print(data['y_op'])
# print(data['y_desc'])
# print(data['y_rotFrame'])
