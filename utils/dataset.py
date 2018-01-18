# Usual tools
import numpy  as np
import pandas as pd

# Combinatorics tools
import itertools

# Root / pandas related tools
from   root_numpy import root2array
import pandautils as pup


class Jet(object):
    """
    Jet Class
    """
    def __init__(self, pt = 0 , eta = 0, phi = 0, m = 0):
        super(Jet, self).__init__()
        
        self._pt = pt
        self._eta = eta
        self._phi = phi
        self._m = m
        
    def __str__(self):
        return 'jet(pT,eta,phi,m)=({:.2f},{:.2f},{:.2f},{:.2f})'.format(self._pt, self._eta, self._phi, self._m)
        
    @property
    def pt(self):
        return self._pt
        
    @property
    def eta(self):
        return self._eta

    @property
    def phi(self):
        return self._phi
        
    @property
    def m(self):
        return self._m
    
    
class Muon(object):
    """
    Muon Class
    """
    def __init__(self, pt = 0 , eta = 0, phi = 0):
        super(Muon, self).__init__()
        self._pt = pt
        self._eta = eta
        self._phi = phi
        
    def __str__(self):
        return 'mu(pT,eta,phi)=({:.2f},{:.2f},{:.2f})'.format(self._pt, self._eta, self._phi)
        
    @property
    def pt(self):
        return self._pt
        
    @property
    def eta(self):
        return self._eta

    @property
    def phi(self):
        return self._phi

    
class Electron(object):
    """
    Electron Class
    """
    def __init__(self, pt = 0 , eta = 0, phi = 0):
        super(Electron, self).__init__()
        self._pt = pt
        self._eta = eta
        self._phi = phi
        
    def __str__(self):
        return 'el(pT,eta,phi)=({:.2f},{:.2f},{:.2f})'.format(self._pt, self._eta, self._phi)
        
    @property
    def pt(self):
        return self._pt
        
    @property
    def eta(self):
        return self._eta

    @property
    def phi(self):
        return self._phi



def add_objectCollection(obj,kinem_name,collection_name,df):
    obj_array = [ obj(pt, eta, phi) for (pt, eta, phi) in zip(
        pup.flatten(df[kinem_name[0]]),
        pup.flatten(df[kinem_name[1]]),
        pup.flatten(df[kinem_name[2]])
    )]
    objects = np.array( pup.match_shape( np.array(obj_array) , df[kinem_name[0]] ) )
    df[collection_name] = objects
    return


def GetdR(j1,j2):
    return np.sqrt( np.square(j1.eta-j2.eta)+ np.square(j1.phi-j2.phi) )


def GetmindR(coll1,coll2):
    if (np.array_equal(coll1,coll2)): 
        all_pairs = list(itertools.combinations(coll1,2))
    else:
        all_pairs = list( itertools.product(coll1,coll2) )
    if(len(all_pairs)>0):
        return np.min([GetdR(o1,o2) for (o1,o2) in all_pairs] )
    else:
        return -1


def compute_manydR(c):
    dRmj=GetmindR(c['muonCollection']     , c['jetCollection'])
    dRjj=GetmindR(c['jetCollection']      , c['jetCollection'])
    dRej=GetmindR(c['electronCollection'] , c['jetCollection'])
    dRee=GetmindR(c['electronCollection'] , c['electronCollection'])
    return dRmj,dRjj,dRej,dRee


def get_data():
    data = pd.DataFrame( root2array('../VectorNtuple_4topSM.root', 'nominal_Loose').view(np.recarray) )
    add_objectCollection(Jet     ,['jet_pt','jet_eta','jet_phi'], 'jetCollection'     , data)
    add_objectCollection(Muon    ,['mu_pt' ,'mu_eta' ,'mu_phi' ], 'muonCollection'    , data)
    add_objectCollection(Electron,['el_pt' ,'el_eta' ,'el_phi' ], 'electronCollection', data)
    return data
