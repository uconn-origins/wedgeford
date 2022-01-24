import os
import shutil
import fnmatch
import sys

def run_RT(output_filename1= 'gas_model.out',curr_name='model',LUV=1.0,LX=1e30, n = 15, p = 300000, r = 200,
          albedo_dir = 'albedo_example'):

    #define values
    s = ' '
    #curr_name_dat = curr_name+".dat" #param filename w/ .dat ext.
    curr_name_uv = curr_name+"_uv"
    curr_name_xray = curr_name+"_xr"

    ####################################################
    #runs both rad transfer UV and rad transfer XRay using terminal commands
    #UV rad transfer
    os.system('./runrt2d_uv_transfer.py -o'+s+curr_name_uv+s+'-i'+s+output_filename1+s+'-n' + s + str(n) + s + '-p' + s + str(p) + s + '-r' + s + str(r) + s ' -E 1,2,3,4,5,7,9 -k ../models/grain_properties/'+ albedo_dir + '/')
    #sys.exit()
    os.system('./runrt2d_xray_transfer.py -o'+s+curr_name_xray+s+'-i'+s+output_filename1+s+'-n ' + s + str(n) +     s + '-p' + s + str(p) + s + '-r' + s + str(r))

    #####################################################
    #runs the CompileUV and CompileXRay codes

    #CompileUV
    os.system("./CompileUVMC.py -n"+s+curr_name_uv+s+"-L" + s +  str(LUV))

    os.system("./CompileXraysMC.py -n "+s+curr_name_xray+s+"-L" + s + str(LX))