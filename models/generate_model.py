from pylab import *
import numpy as np

class grid:
    def __init__(self):
        self.edges = {}
        self.cells = {}
        self.ncell = {}
        self.dx = {}

    def set_grid(self,dim = 'r', xmin=0.1, xmax=500, nx= 100, log= True,input=None):
        valid_dim = {'r','theta','phi'}
        if dim not in valid_dim:
            raise ValueError('valid inputs for dim are r, theta, or phi')
        if input is None:
            if log == True:
                xmin,xmax = log10(max(xmin,0.01)),log10(xmax)
                self.edges[dim] = np.logspace(xmin,xmax,nx+1)
            else:
                self.edges[dim] = np.logspace(xmin,xmax,nx+1)
        else:
            self.edges[dim] = np.sort(np.flatten(np.unique(input)))
        self.ncell[dim] = len(self.edges[dim]-1)
        self.dx[dim] = np.diff(self.edges[dim])
        self.cells[dim] = (self.edges[dim][1:] + self.edges[dim][:-1])*0.5

    def refine_grid(self,dim='r',range_refine=[10,100],factor=4):
        valid_dim = {'r','theta','phi'}
        if dim not in valid_dim:
            raise ValueError('valid inputs for dim are r, theta, or phi')
        valid_factor = int
        if type(factor) != valid_factor or factor < 1:
            raise ValueError('factor should be an integer greater than 1')
        values = self.edges[dim]
        refine_zone = np.arange(0,self.ncell[dim],1)[(values >= range_refine[0]) & (values < range_refine[1])]
        print(refine_zone)
        for i in refine_zone:
            values[i] = np.linspace(values[i],values[i+1],factor)
        values = np.sort(np.unique(sum(values)))

    def make_grid(self,dims=['r','theta']):
        return np.meshgrid(self.cells['r'],self.cells['theta'])


class infall:
    def __init__(self,grid_obj):
        if isinstance(grid_obj,grid) != True:
            raise TypeError('infall needs a grid class object as input')
        elif ['r','theta'] not in grid_obj.cells.keys():
            raise ValueError('grid does not have required dimesions r/theta')
        self.R_C = None
        self.R_IN = None
        self.M_IN = None
        self.grid = grid_obj
        self.r,self.theta = self.grid.make_grid()
        self.u = np.cos(self.theta)

    def solve_u0(self,r=None, u=None,r_c=None):
        if r is None:
            r = self.r
        if u is None:
            u = self.u
        if r_c is None:
            r_c = self.R_C
        r_ = r/r_c
        a = (np.sqrt(729. * u**2 * r_**2 + 4 * (3 * r_ - 3)**3) + 27. * u * r_)**(1./3.)
        b = 2**(1/3.)
        self.u0 = -b * (-1. + r_) / a + a / 3. / b

    def calc_density(self,u0=None,r=None,u=None,r_c = None,m_in=None):
        if u0 is None:
            u0 = self.u0
        if self.u0 is None:
            self.solve_u0(r,u,r_c)
        if r is None:
            r = self.r
        if u is None:
            u = self.u
        if r_c is None:
            r_c = self.R_C
        r_ = r/r_c


g=grid()
g.set_grid('r')
g.set_grid('theta',xmin=0,xmax=pi/2,log=False)
