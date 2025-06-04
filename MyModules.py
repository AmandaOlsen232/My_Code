from matplotlib import rcParams
from cycler import cycler
import numpy as np

class MyPlot:
    def __init__(self):
        # Define your preferred settings
        self.settings = {
            'font.size': 20,
            'font.family': 'Times New Roman',
            'axes.labelsize': 20,
            'axes.titlesize': 16,
            'legend.fontsize': 14,
            'xtick.labelsize': 18,
            'ytick.labelsize': 18,
            
            'lines.linewidth': 2.5,
            'axes.linewidth': 2,
            'axes.prop_cycle': cycler(color=['black']),
            
            'figure.figsize': (8, 6),
            'savefig.dpi': 300,
            'figure.dpi': 100,
            
            'axes.grid': True,
            'grid.alpha': 1.0,
            'grid.color': 'gray',
            'grid.linestyle': '--',
            'grid.linewidth': 1
        }
        
        rcParams.update(self.settings)
        
    def update(self):
        rcParams.update(self.settings)

def jacobian(f, x, h=1e-3):
    '''

    Parameters
    ----------
    f : function
        name of a function that takes in the state vector, x
        the output has dimension, m.
    x : nx1 array
        n dimensional state vector.

    Returns
    -------
    J : array
        returns mxn Jacobian matrix.

    '''
    
    dim_x = len(x)
    dim_fx = len(f(x))
    J = np.zeros((dim_fx, dim_x)) #initialize Jacobian to the appropriate dimension
    
    #loop through the functions
    for j in range(dim_fx):
        #loop through the input variables
        for i in range(dim_x):
            eps = np.zeros(dim_x).reshape(dim_x,1)
            eps[i][0] = h
            
            fp = f(x+eps)
            fm = f(x-eps)
            J[j][i] = (fp[j][0] - fm[j][0])/(2*h) #central differencing
    return J
        
def multivariable_newtons_method(f, x):
    #currently, this function only allows for a symmetric Jacobian (so len(x)=len(f(x)))
    #in the future, I could try and apply a pseudo inverse or something
    #I also need to update the convergence criteria
    '''
    
    Parameters
    ----------
    f : function
        name of a function that takes in the state vector, x
        the output has dimension, m.
    x : nx1 array
        n dimensional state vector.

    Returns
    -------
    x : array
        resulting state vector.
    r : array
        resulting residuals.
    '''
    #do multivariable Newton's method
    for _ in range(4):
        J = jacobian(f, x)
    
        x = x - np.linalg.inv(J) @ f(x)

    r = f(x)
    return x, r

