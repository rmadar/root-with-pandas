######################################
#
# For ~20k Evts:
#
#Total delta  : 10.24s
# --> loading : 1.51s
# --> histo   : not-defined
# --> plotting: 8.73s
#
#####################################

# Timing management
from timeit import default_timer

# Usual numeric & data tools
import itertools as itr
import numpy     as np
import pandas    as pd

# Plotting tools
import matplotlib        as mpl
import matplotlib.pyplot as plt
mpl.rcParams['legend.frameon' ] = False
mpl.rcParams['legend.fontsize'] = 'xx-large'
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['axes.titlesize' ] = 18
mpl.rcParams['axes.labelsize' ] = 18
mpl.rcParams['lines.linewidth'] = 2.5

# ROOT Dataset manipulation
from   root_numpy import root2array
import pandautils as pup
import dataset    as d


#==========================================
#=========    Some functions   ============
#==========================================

def add_objectCollection(obj,kinem_name,collection_name,df):
    """
    This function add a row with an object obj (Jet, Muon or Electron)
    to the dataframe df. Since the current objects need pt,eta,phi
    (stored as vector in the current dataframe), corresponding string
    need to be passed.
    """
    obj_array = [ obj(pt, eta, phi) for (pt, eta, phi) in zip(
        pup.flatten(df[kinem_name[0]]),
        pup.flatten(df[kinem_name[1]]),
        pup.flatten(df[kinem_name[2]])
    )]
    objects = np.array( pup.match_shape( np.array(obj_array) , df[kinem_name[0]] ) )
    df[collection_name] = objects
    return


def GetdR(o1,o2):
    """
    Compute euclidien distance in (eta,phi) space between the
    two objects o1 and o2.
    """
    return np.hypot( o1.eta-o2.eta, o1.phi-o2.phi )


def GetmindR(coll1,coll2):
    """
    Return the minimum dR among all computed on all possible pairs (o1,o2)
    where o1 and o2 are in the collection coll1 and coll2 respectively.
    """
    if (np.array_equal(coll1,coll2)): all_pairs = list( itr.combinations(coll1,2) )
    else:                             all_pairs = list( itr.product(coll1,coll2)  )
    if(len(all_pairs)>0): return np.min([GetdR(o1,o2) for (o1,o2) in all_pairs] )
    else:                 return -1
#===========================================




#======================
# 1. Get and shape data
#======================
t_0 = default_timer()

needed_branches = ['jet_pt','jet_eta','jet_phi','mu_pt' ,'mu_eta' ,'mu_phi','el_pt' ,'el_eta' ,'el_phi' ]
data_array=[]
Nfiles=5
for i in range(0,Nfiles):
    thisdata = pd.DataFrame( root2array('../VectorNtuple_4topSM.root', 'nominal_Loose', branches=needed_branches).view(np.recarray) )
    add_objectCollection(d.Jet     ,['jet_pt','jet_eta','jet_phi'], 'jetCollection'     , thisdata)
    add_objectCollection(d.Muon    ,['mu_pt' ,'mu_eta' ,'mu_phi' ], 'muonCollection'    , thisdata)
    add_objectCollection(d.Electron,['el_pt' ,'el_eta' ,'el_phi' ], 'electronCollection', thisdata)
    data_array.append(thisdata)
data=pd.concat(data_array)
data.info()
t_load = default_timer()


#================
# 2. Compute dRs
#================
dRmj = data.apply( lambda r: GetmindR(r['muonCollection']    , r['jetCollection']     ), axis=1)
dRjj = data.apply( lambda r: GetmindR(r['jetCollection']     , r['jetCollection']     ), axis=1)
dRej = data.apply( lambda r: GetmindR(r['electronCollection'], r['jetCollection']     ), axis=1)
dRee = data.apply( lambda r: GetmindR(r['electronCollection'], r['electronCollection']), axis=1)
t_dR = default_timer()

# ===========
# 3. Plot it
# ===========
plt.figure(figsize=(10,7))
hist_args = {'bins':100, 'alpha':0.8, 'density':True, 'histtype':'step', 'linewidth':3}
ax=plt.hist(dRmj[dRmj>0], label='$dR(\mu,j)$', **hist_args)
ax=plt.hist(dRej[dRej>0], label='$dR(e,j)$'  , **hist_args)
ax=plt.hist(dRjj[dRjj>0], label='$dR(j,j)$'  , **hist_args)
ax=plt.hist(dRee[dRee>0], label='$dR(e,e)$'  , **hist_args)
plt.xlabel('min$_{\{i,j\}}$$\{dR(obj_i, obj_j)\}$')
plt.ylabel('PDF')
plt.yscale('log', nonposy='clip')
plt.legend()
plt.tight_layout()
plt.savefig('dR.png')
t_plot = default_timer()


# ===============
# 4. Print result
# ===============
print( '\nTotal Number of events     : {:.0f}'.format(len(data)) )
print( 'Total delta                : {:.2f}s'  .format(t_plot-t_0) )
print( ' --> loading               : {:.2f}s'  .format(t_load-t_0) )
print( ' --> dR computation        : {:.2f}s'  .format(t_dR-t_load) )
print( ' --> plotting and selection: {:.2f}s'  .format(t_plot-t_dR) )

