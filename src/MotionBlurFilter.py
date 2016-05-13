
# Written by Chris Calderon 2015 (Chris.Calderon@UrsaAnalytics.com)
#
#
# Copyright 2015 Ursa Analytics, Inc.

   # Licensed under the Apache License, Version 2.0 (the "License");
   # you may not use this file except in compliance with the License.
   # You may obtain a copy of the License at

   #     http://www.apache.org/licenses/LICENSE-2.0
   
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

__license__ = "Apache License, Version 2.0"
__author__  = "Chris Calderon, Ursa Analytics, Inc. [www.UrsaAnalytics.com]"
__status__  = "Development"

import numpy as np
import scipy.special as spspecial




class ModifiedKalmanFilter1DwithCrossCorr(object):
    def __init__(self,tsData,dt,StaticErrorEstSeq=None): 
        """
        setup class for "Motion Blur filter" for 1D uniformly sampled case where state evolves in continuous time 
        SDE driven by standard Brownian motion and discrete measurements are "blurred"  
        (blurring occurs both from "dynamic" motion blur and "static" point spread function fitting errors).   
        filtering code ignores computational issues inherent to multivariate models.  
        also assume stationary distribution for initial cov vs. information filter.

        Note: although code can handle estimation cases with \kappa near zero, recommended to use MA1 code provided if \kappa is within
        parameter uncertainty of zero since likelihood expansion used for limiting case (analytic limit exists, this numerical implementation 
        just switches to taylor likelihood proxy to avoid potential numerical overflow (this introduces some slight bias in estimates obtained for very small kappa).   
        if true  \kappa < 0 in DGP, model isn't "confined" and current code should not be used for this "unstable" case).   improved algorithm could be made to "take limits" vs. taylor proxy;
        i.e., commercial quality code can be made to avoid taylor likelihood expansion, but this code is useful for illustrating basic ideas. 

        class setup to allow "traditional" kalman filter by redefining some methods and using some aux variables (modifications for 1D case illustrated below) 

        input:
        tsData:  np.array of time series data
        dt:      float giving time spacing between observations.  straightforward to adapt to nonuniform time spacing (just need to compute filter parameters and/or input time changing R_k sequence)
        StaticErrorEstSeq: [optional].  list or np.array of same length as tsData with STD DEV est of static errors.  allows fusing in estimations (default zero)
                           if this parameter is input, the addition of the estimated \hat{sigma}^{loc} to the input stream gives a refined estimate of the empirical "static noise" 

        output: 
        likelihood: log likelihood of innovations
        xfilt:      filtered state for given parameters  
        pit:        PIT associated with likelihood function evaluated at given parameters
        Shist:      history of innovation covariance evaluated at given parameters.  
        """
        
        tsData = np.array(tsData) #cast in case py list passed in (assign column vector shape later)
        self.__T=max(tsData.shape) #compute length of time series
        self.__y=np.reshape(tsData,(self.__T,1),order='F')#reshape data to column vector (scalar time series assumed)
        self.__dt=dt #store internally to simplify interface to 3rd party opt functions

        if StaticErrorEstSeq is None:
            self.__Rbase = [0]*self.__T #initialize with zeros.
            
        else:
            if len(StaticErrorEstSeq) == self.__T:
                self.__Rbase =[loci for loci in StaticErrorEstSeq]
                
            else:
                print 'WARNING:'*5
                print 'Input StaticErrorEstSeq has different length the observation sequence.'
                print 'Using first entry in list for entire duration for computations.'
                print 'WARNING:'*5
                self.__Rbase = [StaticErrorEstSeq[0]]*self.__T #quantity always squared, so need to ensure input is "std. loc estimate" and not variance  

    def KFfilterOU1d(self,pars,evalInnovStats='True',P10=None):
        #run naive simple filter where scalar is observed.  
        #note scalar case doesn't need to worry about symm pos. def. mats and other computational issues commonly encountered in multivariate case  
        #***** make computations easy to follow from a "textbook" classic KF implementation;  
        #***** a module executing the "standard" Kalman filter where correlation between process and measurement noise is zero \forall t is provided to illustrate connection
        #***** between motion blur filter and standard Kalman filter

        x10 = np.mean(self.__y[0:5]) #initialize with mean of first few observations (ignore smoothing and information filter approaches for simple illustrative code)
        delta = self.__dt
        kappa = pars[0]
        sig   = pars[1]

      
        
        F,Q,H,R,C,HF,A,HA = self.PttUpdateOU1dPars(pars,delta)
        
        if P10 is None:  # assume stationary dist instead of info filter for initial uncertainty if no user supplied uncertainty provided (zero uncertainty in state position  inference poses a serious technical problem in MLE anyway...)
            if np.abs(kappa)>1E-5:
                P10 = sig**2/kappa/2. #standard stationary confined diffusion covariance
            else:
                P10 = np.abs(delta*sig**2 - delta**2*kappa*sig**2) #if small kappa is passed in, use a small t expansion (close to diffusion) to permit 
                #directed and pure diffusion models without divide by zero issues (taking kappa -> 0 gives exact directed or pure diffusion likelihood, 
                #but don't trust computational routine to correctly  "infer" limits)

        P_Innov = P10
        x00 =  x10 #keep explicit copy of filtered estimate from previous iteration for special blurred filter 
     
 
        xfilt =[] #track the filtered estimates
        pit =[]  #return the pit random variables
        Shist = []
        loglikelihood = 0.
        Rconstant = R  #store a copy of the value dictated by the MLE parameter (adjust R_t each cycle of the computation)
        

        for idx, y in enumerate(self.__y):
 
            R = Rconstant + ( self.__Rbase[idx] + pars[2] )**2 #use the input static error sequence (squared) to give time dependent R.
            

            if evalInnovStats:
                #compute innovation statistics
                Sinv = 1./(H*P_Innov*H+R) #use variable P_Innov to stress fundamental difference in innovation computation 
                #between motion blur and classic kalman filter
                Sinv = max((Sinv,np.finfo(float).resolution)) #time varying feature or roundoff can introduce negative innovation variance.  only permit values above machine precision
                                                              #if MLE has Sinv < 0, GoF tests will readily identify this issue. 
                z = (y-HF*x00-HA)*np.sqrt(Sinv) #HF and HA are other nuances of blurred filter formulation used 
                piti = spspecial.erfc(z/np.sqrt(2.))*.5 #compute CDF of normal
                loglikelihood += 1/2.*(np.log(Sinv)) + -z*z/2. -1/2.*np.log((2*np.pi))
                pit.append(piti)
                Shist.append(1./Sinv)
               
            
            #compute gain and then fuse in information from current measurement 
            K = self.computeGain(P_Innov,C,H,R,F) #different gain computation for "classic" KF and MBF
            # K = self.computeGain(P10,C,H,R,F)
            x11 = x10 + K*(y-HF*x00-HA) #HF and HA are nuances of blurred filter formulation used 
            x00 = x11 #keep explicit copy of filtered estimate from current iteration for special blurred filter likelihood eval
            xfilt.append(x11) #store filtered estimate
            
            #update/forecast state for simple mean zero OU model
            x10=F*x11 + A
            
            P00 = P10 -  K*(H*P_Innov*H+R)*K 
            P10 = F*P00*F  +  Q

            P_Innov = self.P_Innov(P10,P00) #nuance of motion blur filter
 
            

        
        xfilt = np.array(xfilt)
        loglikelihood = loglikelihood/self.__T #return empirical time average of loglikelihood

        return loglikelihood,xfilt,pit,Shist #return all stats (create wrapper to make more efficient feval calls)
        

    def evalCostFunc(self,pars):
        feval = self.KFfilterOU1d(pars,evalInnovStats='True')
        negloglike = -feval[0][0]
        return negloglike #return negative loglikehood for minimization routines (also set flags to make computation more efficient) 

    def computeGain(self,P10,C,H,R,F):
            K   = (C+F*P10*H)/(H*P10*H+R) #blur form for updating covariance of filtered state.
            return K

    def P_Innov(self,P10,P00): #simple switch function permitting both the classic and motion blur filter with one code base
            return P00 #blur form for updating covariance of filtered state.
            

    def PttUpdateOU1dPars(self,pars,delta):
        #par is assumed to contain (kappa,sigma,stdloc)
        kappa = pars[0]
        sigma = pars[1]
        R     =  0 #in time varying code, permit negative parameters reducing input variance and assign 
                   #localization contribution to net measurement noise  in main code (here just assign "blur" contribution)
        

        if len(pars)>3:
            alpha = pars[3]
        else:
            alpha = 0

        F     = np.exp(-kappa*delta) #standard res
        
        #In order to avoid numerical problems with machine zero kappa ("pure directed diffusion"), 
        #use taylor proxies if kappa is around machine single precision.  keep standard KF result simple and unable
        #to handle this special case.
        if np.abs(kappa)>1E-5:
            Q     = (sigma**2/2./kappa)*(1.-np.exp(-2.*kappa*delta))
            Qblur = 2*delta/kappa - 3/kappa**2 + 4*np.exp(-delta*kappa)/kappa**2 - np.exp(-2*delta*kappa)/kappa**2
            Qblur = Qblur*(sigma**2/2./kappa)/(delta**2)
            H     = (1. - np.exp(-delta*kappa) )/kappa/delta
            #compute the exact cross correlation term of time integrated OU vs. discretely sampled state (S in notation of Anderson and Moore...I prefer using S for innovation covariance)
            C     = (1./kappa - 2.*np.exp(-delta*kappa)/kappa + np.exp(-2.*delta*kappa)/kappa)*sigma**2/kappa/2./delta

            fp = alpha/kappa 
            
            A  = (1-F)*fp #form assumes kappa>0 implies stable linear system
        #compute Integral((1-exp(-kappa*(s)))*alpha/delta/kappa,(s,0,delta)) [form also assumes kappa>0 implies stable linear system]
            HA =  fp - fp/(delta*kappa) + fp*np.exp(-delta*kappa)/(delta*kappa)
            
        else: #note: expansion results not thoroughly tested.  if kappa is smaller than Cramer Rao asymptotic bound, 
              #recommended to use MA(1) model code provided with constant "velocity" adjustment for exact likelihood for statistical inference (point estimates expected to be okay, but GoF and other statistics requiring higher likelihood accuracy questionable) 
            A  = alpha*(delta - kappa*delta**2/2.) # + O(delta^3) + O(kappa^2)
            HA = alpha*delta/2 - alpha*delta**2*kappa/6 + alpha*delta**3*kappa**2/24 # + O(delta^4) + O(kappa^3)
            Q  = delta*sigma**2 - delta**2*kappa*sigma**2
            Qblur = 2*delta**1/3. - delta**2*kappa/2. + 7*delta**3*kappa**2/30.
            #compute expansion (in kappa) of cross correlation term of time integrated OU vs. discretely sampled state (S in notation of Anderson and Moore)
            C  = delta*sigma**2/2 - delta**2*kappa*sigma**2/2 + 7*delta**3*kappa**2*sigma**2/24. 
            H  = 1. 
        
        R     += Qblur #add blur contribution to effective measurement noise 
        HF    = H #special case for blur model
        
        return F,Q,H,R,C,HF,A,HA 


class ClassicKalmanFilter(ModifiedKalmanFilter1DwithCrossCorr):
    """
    generates parameters for using the "blur" version of the 1D KF filter with the "classic Kalman filter" where there is no
    statistical dependence / correlation between process and measurement noise.

    for I/O and notes, see parent class.  the methods redefined here show how to introduce variables and redefine quantities
    to implement the "classic" KF.
    """
    def __init__(self,tsData,dt,StaticErrorEstSeq=None):
        super(ClassicKalmanFilter, self).__init__(tsData,dt,StaticErrorEstSeq)
        
    def computeGain(self,P10,C,H,R,F):
        K   = (P10*H)/(H*P10*H+R) #gain form required for using classic KF within "motion blur filter" formulation. 
        #note: C=0 required for "standard" classic KF (enforced in code)
        return K

    def P_Innov(self,P10,P00): #simple switch function permitting both the classic and motion blur filter with one code base
            return P10 #blur form for updating covariance of filtered state.

    def PttUpdateOU1dPars(self,pars,delta):
        #par is assumed to contain (kappa,sigma,stdloc)
        kappa = pars[0]
        sig   = pars[1]
    
        R     = 0 #in time varying code, permit negative parameters reducing input variance and assign 
                   #localization contribution to net measurement noise  in main code (here just assign "blur" contribution)
        if len(pars)>3:
            alpha = pars[3]
        else:
            alpha = 0
        #Keep expression below simple.  just note numerical issues may arise if kappa near machine zero is attempted (practically not too relevant since MA1 case and KF should give numerically identical/similar results)
        F     = np.exp(-kappa*delta) 
        Q     = (sig**2/2./kappa)*(1.-np.exp(-2.*kappa*delta))
        H     = 1.
        HF    = H*F #this generates a standard KF by a means that fits into the MBF framework.  multiplication by F maps ri|i,Pi|i to ri+1|i,Pi+1|i and then mult by H gives observation stats
        fp    = alpha/kappa 
        A     = (1-F)*fp #assumes kappa>0 implies stable linear system
        HA    = H*A # similar to HF above (for classic KF, simply multiplies H by A of state space model)
        C     = 0.
        
        
        
        return F,Q,H,R,C,HF,A,HA


def simTimeIntegratedOUdeltaICvel(pars,delta,Tsub,Nsim,T,xIC=0):
    """
    inputs:

    par   : assumed to contain (kappa,sigma,stdloc,vel)
    delta : time interval between observations
    Tsub  : integer specifying granularity of division between delta (used for simple numerical time integration)
    Nsim  : number of paths to simulate
    T     : length of each path
    xIC   : dirac point mass of ensemble.  code can be modified to handle distribution of ICs., though not relevant to this study.

    outputs:

    xraw  : TxN nparray of (exact) discretely sampled OU process 
    xblur : "     "   simple time quadrature of xraw with accuracy dictated by Tsub (avoid using analytical expression in order to show convergence to numerical results)

    [both outputs check out when compared to analytical expressions]    
    """
    

    kappa = pars[0]
    sig   = pars[1]
    R     = pars[2]**2
    if len(pars)>3:  #allow for varying length inputs (assume constant vel is zero typically)
        vel   = pars[3]
    else:
        vel   = 0.
 
    dt    = float(delta)/float(Tsub)
    F     = np.exp(-kappa*dt)

    if np.abs(kappa)>1E-5:    #this code chunk allows an SPT legacy model "pure directed diffusion" (if kappa^2 is near double precision zero, no computational difference)
        fp = vel/kappa 
        A  = (1-F)*fp #assumes kappa>0 implies stable linear system
        sqrtQ = np.sqrt((sig**2/2./kappa)*(1.-np.exp(-2.*kappa*dt)))
    else:
        A  = vel*(dt - kappa*dt**2/2.) # + O(delta^3) + O(kappa^2)
        sqrtQ = np.sqrt(dt*sig**2 - dt**2*kappa*sig**2) # + O(delta^3)++ O(kappa^2)

    

    
    x0    = np.ones(Nsim)*xIC
    
    xraw=[]
    xblur=[]
    xrawsampled=[]

    for i in range(T):
        xloc=[]
        for j in range(Tsub):
            noise = sqrtQ*np.random.randn(Nsim)
            x0=F*x0+noise+A
            xloc.append(x0)
            xraw.append(x0)

        xrawsampled.append(x0)
        xblur.append(np.mean(xloc,axis=0))

    xraw = np.array(xraw)
    xrawsampled = np.array(xrawsampled)
    xblur= np.array(xblur)
    

    return xraw,xrawsampled,xblur 








