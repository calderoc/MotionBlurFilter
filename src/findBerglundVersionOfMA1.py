#!/usr/bin/env python

#collection of routines for estimating diffusion + noise model from different of uniformly sampled observations.
#modification to non-uniform sampling straight-forward (just modify diagonal...this makes eigen-values harder to compute,
#but efficient sparse solvers can be leveraged)

#this versions allows for arbitrary correlation in the diagonal term by using a nonstandard parameterization of 
#the MA1.  nonstandard parameterization useful when blur dominates in magnitude (i.e., larger than static and diffusive noise).


import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as sla
import scipy.sparse as sparse
import timeit


class CostFuncMA1Diff(object):
    def __init__(self,tsData,dt,blurCoef=1./6.): 
        """
        class for setting up costfunc of MA1 of differenced measurements in diffusion plus noise model
        dx_t= v dt+ sqrt(sqrtD*2)dBt
        y_ti = x_ti + \epsilon_i*sqrtR
        pars:=(sqrtD,sqrtR,v)
        dt:= scalar giving time sampling (can modify to vector dt if need be...see comments for code changes required)
        blurCoef [optional].  parameter corresponding to "R" in Berglund's 2010 PRE paper (defaults to uniform continuous illumination value)
        """
        

        self.__T=max(tsData.shape)-1 #compute length of differenced time series
        tsData=np.reshape(tsData,(self.__T+1,1),order='F')#reshape data to column vector (scalar time series assumed)
        self.__dy=np.diff(tsData,axis=0) #compute differenced time series (do not manipulate this, modify local copies via -pars[2]*dt)
        
        ii=np.arange(self.__T) #construct and store coo sparse indices for MA1
        self.__rowi=np.concatenate((ii,ii[:-1],ii[:-1]+1))
        self.__coli=np.concatenate((ii,ii[:-1]+1,ii[:-1]))
        self.__dt=dt #store internally to simplify interface to 3rd party opt functions
        self.__blurCoef=blurCoef

    def evalCostFuncVel(self,pars):
        """
        interface method to par opt routine
        MA1 of differenced measurements in diffusion plus noise model
        dx_t= v dt+ sqrt(sqrtD*2)dBt
        y_ti = x_ti + \epsilon_i*sqrtR
        pars:=(sqrtD,sqrtR,v)
        dt:= scalar giving time sampling (can modify to vector dt if need be...see comments for code changes required)
        """
        
        #use Berglund parameterization (allows for negative or positive correlation in MA1 covariance)
        Reff   = pars[1]**2 - 2.*pars[0]**2*self.__dt*self.__blurCoef #use equation 7 from Berglund for R = 1/6. (permits negative evals) 
        Deff   = pars[0]**2*self.__dt

        #compute MA1 covariance matrix given effective parameters above
        S = self.constructMA1(Deff,Reff)


        valsE  = np.array([1.]*self.__T)*2.*Deff #formulation assume covariance of form: 2*Deff*Id + Reff*SecondDiffMat
        tmp    = np.arange(self.__T)+1
        valsE2 = (4.*Reff)*  ((np.sin(tmp*np.pi/2./(self.__T+1)))**2) #pure dirichlet boundary conditions
        vals   = valsE + valsE2 #

        loglikelihood=sum(np.log(vals))/2.
        
        #compute quadratic form contribution to log likelihood
        dy=self.__dy-pars[2]*self.__dt #
        #execute solve required to compute quadratic form
        tmp = self.sparseSolveWrapper(S,dy)
        quadForm = np.dot(dy.T,tmp)

        loglikelihood+=quadForm/2. #TODO:  replace this line by adding in quadratic form.
        #note negative of (unormalized) loglikelihood  computed above 


        return loglikelihood


    
    def constructMA1(self,Deff,Reff):
        """
        precompute the coo sparse matrix indices of a tri-banded MA1  matrix (stored in rowi, coli) and return sparse coo mat
        pars:=(sqrtD,sqrtR,v)
        dt:= scalar giving time sampling (can modify to vector dt if need be)
        """
        #form sparse MA1 matrix
        R=Reff
        mainDiag=(2*R+2*Deff)*(np.array([1.]*self.__T)) #expression "np.array([1.]*N" like matlab ones(N,1) (with row/column left open)
        band=-R*(np.array([1.]*(self.__T-1)))
    
        svals=np.concatenate((mainDiag,band,band))
        svals=np.array([float(i) for i in svals]) #crude approach to computing a array with shape (T,) vs (T,1).  difference required for sparse
        S=sparse.coo_matrix((svals,(self.__rowi,self.__coli)),shape=(self.__T,self.__T))
        return S
    
    def sparseSolveWrapper(self,S,RHS):

        supalu = sla.splu(S.tocsc())
        tmp = supalu.solve(RHS.reshape(-1))
         
        return tmp


def runOPT(sig1=.4,sig2=.4,sqrtR=35/100.,dt=10/1000.,N=10,T=50):
    
    ts=np.arange(N)

    dW,W=simSDE.simDB(T=T,N=N,dt=dt,sig=sig1*np.sqrt(2))
    dW2,W2=simSDE.simDB(T=T,N=N,dt=dt,sig=sig2*np.sqrt(2))
    

    W2+=np.reshape(W[:,-1],(N,1)) #create a smooth transition by adding terminal value of 
    print W.shape
    W=np.hstack((W,W2))
    print W.shape

    Y=W+np.random.randn(W.shape[0],W.shape[1])*sqrtR
    fracsplit=.5 #adjust this parameter to reflect mix of sig1 and sig2 in sampled data
    sigEff = np.sqrt(fracsplit*2*sig1**2+fracsplit*2*sig2**2) 
    sigEff = sigEff/2. #make sigEff^2= D_Eff 
    Xtrue= np.array([sigEff,sqrtR,0])
    
    #iterate over paths and carry out optimization
    resH=[]
    for i,yi in enumerate(Y):
        print 'Iteration:',i
        costInstance = CostFuncMA1Diff(yi,dt)
        res = spo.minimize(costInstance.evalCostFunc, Xtrue/2., method='nelder-mead',options={'xtol': 1e-5, 'disp': False})
        print res.x
        resH.append(res.x)
    resH=np.asarray(resH)

    print '******* Result Summary ***********************'
    print ''
    print 'True (or Effective) Par:', Xtrue
    print ''
    print 'parameter means,medians, max, min of NxPar history vec:'
    print np.mean(np.abs(resH),axis=0) #takes abs value since optimization was unconstrained (cost function squares sig and sqrtR, so no diff;  physically both pars must be >0)
    print np.median(np.abs(resH),axis=0)
    print np.max(np.abs(resH),axis=0)
    print np.min(np.abs(resH),axis=0)
    print 'parameter STD of NxPar history vec:'
    print np.std(np.abs(resH),axis=0)

    print 'Ddt/2R:', sig1**2*dt/(2.*sqrtR**2)
    fname = '~/optRes_dt_'+str(int(dt*1000)) + '_sig1_' + str(sig1) +  '_sig2_' + str(sig1) +'.h5'
    print "saving  results to file: ", fname
    wouth5(fname,resH)

def simDB(N=1000,T=20,dt=2./100.,sig=1): #write out an ASCII file to "floc" given a 2D numpy array "DATAMAT"
    """
      simulate batch of increments and paths of pure BM.  
    """
    sdt=np.sqrt(dt)
    dW=np.random.randn(N,T)*sdt*sig;
    W=np.cumsum(dW,axis=1)
    return dW,W

def MSD(ymeas,tlag=None):
    """
    simple routine for computing MSD.  note overlapping intervals permitted (nonoverlapping intervals noisier, but "uncorrelated" under standard assumptions)
    """
    T = len(ymeas)
    MSD = []
    if tlag is None:
        tlag = range(T/10)
    for t in tlag:
        tmp = []

        for ix,yi in enumerate(ymeas[:-t]):
            dy2i = (yi-ymeas[ix+t])**2
            tmp.append(dy2i)

        MSD.append(np.mean(tmp))

    return tlag,MSD




 
