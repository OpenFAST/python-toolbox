''' Class for defining linear model 

Will have an A, B, C, and D matrix for every operating point
u_op, y_op, and x_op

'''

import numpy as np
import scipy as sp
import control as co
import pyFAST.linearization.linearization as lin
import matplotlib.pyplot as plt


class LinearTurbineModel(object):


    def __init__(self,fromMat=False, matDict=[]):

        if fromMat:   # from matlab .mat file m
            print('here')

            # operating points
            # u_ops \in real(n_inps,n_ops)
            u_ops = np.zeros(matDict['u_ops'][0][0].shape)
            for u_op in matDict['u_ops'][0]:
                u_ops = np.concatenate((u_ops,u_op),axis=1)

            self.u_ops = u_ops[:,1:]

            # y_ops \in real(n_outs,n_ops)
            y_ops = np.zeros(matDict['y_ops'][0][0].shape)
            for y_op in matDict['y_ops'][0]:
                y_ops = np.concatenate((y_ops,y_op),axis=1)

            self.y_ops = y_ops[:,1:]

            # x_ops \in real(n_states,n_ops), note this is un-reduced state space model states, with all hydro states
            x_ops = np.zeros(matDict['x_ops'][0][0].shape)
            for x_op in matDict['x_ops'][0]:
                x_ops = np.concatenate((x_ops,x_op),axis=1)

            self.x_ops = x_ops[:,1:]

            # Matrices

            # A \in real(n_states,n_states,n_ops)
            n_states = np.shape(matDict['A'][0][0])[0]
            A_ops   = np.zeros((n_states,n_states,1))
            for A_op in matDict['A'][0]:
                A_ops   = np.concatenate((A_ops,np.expand_dims(A_op,2)),axis=2)

            self.A_ops = A_ops[:,:,1:]

            # B \in real(n_states,n_inputs,n_ops)
            n_states = np.shape(matDict['B'][0][0])[0]
            n_inputs = np.shape(matDict['B'][0][0])[1]
            B_ops   = np.zeros((n_states,n_inputs,1))
            for B_op in matDict['B'][0]:
                B_ops   = np.concatenate((B_ops,np.expand_dims(B_op,2)),axis=2)

            self.B_ops = B_ops[:,:,1:]

            # C \in real(n_outs,n_states,n_ops)
            n_states = np.shape(matDict['C'][0][0])[1]
            n_outs = np.shape(matDict['C'][0][0])[0]
            C_ops   = np.zeros((n_outs,n_states,1))
            for C_op in matDict['C'][0]:
                C_ops   = np.concatenate((C_ops,np.expand_dims(C_op,2)),axis=2)

            self.C_ops = C_ops[:,:,1:]

            # D \in real(n_outs,n_inputs,n_ops)
            n_states = np.shape(matDict['D'][0][0])[1]
            n_outs = np.shape(matDict['D'][0][0])[0]
            D_ops   = np.zeros((n_outs,n_inputs,1))
            for D_op in matDict['D'][0]:
                D_ops   = np.concatenate((D_ops,np.expand_dims(D_op,2)),axis=2)

            self.D_ops = D_ops[:,:,1:]

            # Save wind speed as own array since that's what we'll schedule over to start
            self.u_h     = self.u_ops[0]

            # Input/Output Indices
            self.ind_fast_inps     = matDict['indInps'][0] - 1
            self.ind_fast_outs     = matDict['indOuts'][0] - 1



    def solve(self,tt,u_h,Plot=True):


        # interpolate system using uh_op = mean(u_h)
        # maybe make this its own function?  wait and see if we need it anywhere else...
        uh_op = np.mean(u_h)

        f_A     = sp.interpolate.interp1d(self.u_h,self.A_ops)
        f_B     = sp.interpolate.interp1d(self.u_h,self.B_ops)
        f_C     = sp.interpolate.interp1d(self.u_h,self.C_ops)
        f_D     = sp.interpolate.interp1d(self.u_h,self.D_ops)

        f_u     = sp.interpolate.interp1d(self.u_h,self.u_ops)
        f_y     = sp.interpolate.interp1d(self.u_h,self.y_ops)
        f_x     = sp.interpolate.interp1d(self.u_h,self.x_ops)

        A       = f_A(uh_op)
        B       = f_B(uh_op)
        C       = f_C(uh_op)
        D       = f_D(uh_op)

        u_op    = f_u(uh_op)
        x_op    = f_x(uh_op)
        y_op    = f_y(uh_op)


        # form state space model
        P_op            = co.StateSpace(A,B,C,D)
        P_op.inputs     = ['WindSpeed','BldPitch']
        P_op.outputs    = ['GenSpeed','TwrBsMyt','PltPitch','NacIMU']

        # linearize input (remove uh_op)
        u_lin       = np.zeros((2,len(tt)))
        u_lin[0,:]  = u_h - uh_op

        _,y_lin,xx = co.forced_response(P_op,T=tt,U=u_lin)

        # Add back in operating points
        y   = y_op[self.ind_fast_outs].reshape(-1,1) + y_lin
        u   = u_op[self.ind_fast_inps].reshape(-1,1) + u_lin
        print('here')

        if Plot:
            ax = plt.subplot(411)
            ax.plot(tt,u[0,:]+uh_op)
            ax.set_ylabel('WindSpeed')
            ax.set_xticklabels([])
            ax.grid(True)

            ax = plt.subplot(412)
            ax.plot(tt,y[0,:])
            ax.set_ylabel('GenSpeed')
            ax.set_xticklabels([])
            ax.grid(True)

            ax = plt.subplot(413)
            ax.plot(tt,y[1,:])
            ax.set_ylabel('TwrBsMyt')
            ax.set_xticklabels([])
            ax.grid(True)

            ax = plt.subplot(414)
            ax.plot(tt,y[2,:])
            ax.set_ylabel('PltPitch')
            ax.grid(True)

            plt.show()
        
        print('here')

