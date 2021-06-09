## Combine output of dales to an input file suitable for tenstream solver 
##
## example usage: srun -p short python create_input_cld_fld.py -procx 16 -procy 24 -path /scratch-shared/<..>/<..>/<..>/ -time 18000 -sim 1 -n_proc -1

import time
import netCDF4 as nc
import numpy as np
import multiprocessing as mp
import numba
import argparse

parser = argparse.ArgumentParser(description="Prepare DALES output for the tenstream solver")
parser.add_argument("-path", type=str, default="./", help="path where data is located")
parser.add_argument("-sim",  type=int, default=1, help="simulation number")
parser.add_argument("-procx", type=int, default=1, help="processors in x direction")
parser.add_argument("-procy", type=int, default=1, help="processors in y direction")
parser.add_argument("-time", type=int, default=0, help="time to save ")
parser.add_argument("-n_proc", type=int, default=1, help="number of processors to use (-1=all available)")
args = parser.parse_args()

path = args.path
sim_name  = "{:03d}".format(args.sim)
basename = path+"threedheating.{:03d}.{:03d}.{}.nc"
basename_sfc = path+"surfcross.x{:03d}y{:03d}.{}.nc"
basename_fld = path+"field.{}"
proc_x = args.procx
proc_y = args.procy
t_save = args.time
            
rho_liq = 1000.
Nc_0 = 200e6
sig_g = 1.34
p0 = 1e3
rcp = .286
cp = 1004.
rlv = 2.53e6

@numba.jit(nopython=True)
def calc_re(lwc):
    return 1.e6*( 3.*( 1.e-3* lwc ) /(4. * np.pi*Nc_0*rho_liq) )**(1./3.) * np.exp(np.log(sig_g)**2 )

@numba.jit(nopython=True)
def calc_tabs(thl,lwc,p0p):
   tabs = thl*(p0p**-rcp) + (1e-3 * rlv/cp) * lwc
   return tabs 

@numba.jit(nopython=True)
def calc_plev(p_lay):
   p_lev = np.zeros(len(p_lay)+1)
   p_lev[1:-1] = (p_lay[1:] + p_lay[:-1])/2.
   p_lev[0] = 2*p_lay[0] - p_lev[1]
   p_lev[-1] = 2*p_lay[-1] - p_lev[-2]
   return p_lev

def load_vars(xxx,yyy, itime, lock):
    nc_in = nc.Dataset(basename.format(xxx,yyy,sim_name))
    lwc = nc_in.variables['ql'][itime] * 1e-5 * 1e3 #g/kg
    lwc = np.swapaxes(lwc, 0, 2)
    x0,x1 = xxx*len(x_local), (xxx+1)*len(x_local)
    y0,y1 = yyy*len(y_local), (yyy+1)*len(y_local)
    
    #obtain pressure
    p_lay = np.loadtxt(basename_fld.format(sim_name))[itime*len(z_lay):(itime+1)*len(z_lay),2]
    p_lev = calc_plev(p_lay)
    
    thl = nc_in.variables['thl'][itime].swapaxes(0,2)
    tabs = calc_tabs(thl, lwc, p0/p_lay[np.newaxis,np.newaxis,:])
    
    nc_in.close()
    
    lock.acquire()
    
    ncf = nc.Dataset(path+"cloud_field.nc","a")
    out_lwc  = ncf.variables['lwc']
    out_lre  = ncf.variables['reff']
    out_plev = ncf.variables['plev']
    out_tabs = ncf.variables['tabs']
    out_tsfc = ncf.variables['tsfc']
    
    out_lwc[x0:x1, y0:y1, :] = lwc
    out_lre[x0:x1, y0:y1, :] = calc_re(lwc)
    out_plev[x0:x1, y0:y1, :] = p_lev[np.newaxis,np.newaxis,:]
    out_tabs[x0:x1, y0:y1, :] = tabs
    out_tsfc[x0:x1,y0:y1] = nc.Dataset(basename_sfc.format(xxx,yyy,sim_name)).variables['tskin'][itime].swapaxes(0,1) * (p0/p_lev[0])**-.286

    ncf.close()
    
    lock.release()

    
if __name__ == '__main__':

    ## Read grid
    nc_in = nc.Dataset(basename.format(0,0,sim_name))
    z_lay = nc_in.variables["zt"][:]
    z_lev = nc_in.variables["zm"][:]
    z_lev = np.append(z_lev, z_lev[-1]+(z_lev[1]-z_lev[0]))
    
    y_local = nc_in.variables["yt"][:]
    x_local = nc_in.variables["xt"][:]
    
    dx = x_local[1] - x_local[0]
    dy = y_local[1] - y_local[0]
    
    nx = proc_x * len(x_local)
    ny = proc_y * len(y_local)
    
    times = nc_in.variables['time'][:]
    
    nc_in.close()
    
    ## Create output file
    ncf = nc.Dataset(path+"cloud_field.nc","w")
    ncf.createDimension("z_lay", len(z_lay))
    ncf.createDimension("z_lev", len(z_lev))
    ncf.createDimension("x", nx)
    ncf.createDimension("y", ny)
    
    nc_z_lay = ncf.createVariable("z_lay", "f8", ("z_lay",))
    nc_z_lev = ncf.createVariable("z_lev", "f8", ("z_lev",))
    nc_x = ncf.createVariable("x", "f8", ("x",))
    nc_y = ncf.createVariable("y", "f8", ("y",))
    
    nc_z_lay[:] = z_lay
    nc_z_lev[:] = z_lev
    nc_x[:] = np.arange(0, nx*dx, dx) + dx/2.
    nc_y[:] = np.arange(0, ny*dy, dy) + dy/2.
    
    nc_lwc = ncf.createVariable('lwc', "f8", ("x","y","z_lay"))
    nc_reff = ncf.createVariable('reff', "f8", ("x","y","z_lay"))
    nc_plev = ncf.createVariable('plev', "f8", ("x","y","z_lev"))
    nc_tabs = ncf.createVariable('tabs', "f8", ("x","y","z_lay"))
    nc_tsfc = ncf.createVariable('tsfc', "f8", ("x","y"))
    ncf.close()   

    n_cpu = mp.cpu_count() if args.n_proc==-1 else args.n_proc
    
    itime = abs(times - t_save).argmin()
    print("saving clouds at t={}".format(t_save) if t_save == times[itime] else "saving clouds at nearest time (t={})".format(times[itime]))

    lock = mp.Manager().Lock()
    a = time.time()
    
    with mp.Pool(n_cpu) as pool:
        res = [pool.apply_async(load_vars, args=(xxx, yyy, itime, lock)) for xxx in range(proc_x) for yyy in range(proc_y)]
        tmp = [r.get() for r in res]
    b = time.time()

    print("elapsed timë:",b-a)
