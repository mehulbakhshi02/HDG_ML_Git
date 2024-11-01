import numpy as np
import subprocess
import math
import time
import fileinput
from scipy.stats import loguniform

#from netgen.read_gmsh import ReadGmsh
#from ngsolve import *
start_time = time.time()

def MakeNetgenMesh():
    mesh = ReadGmsh("output")
    mesh = Mesh(mesh)
    Draw(mesh)
    mesh.ngmesh.Save("netgen.vol")
    subprocess.call('cd Scripts; g++ netgen2mmg2d.cpp -o netgen2mmg2d; ./netgen2mmg2d', shell=True)

def DeleteFiles():

    # Commented by KUNAL
    subprocess.call('rm bamg_prev.mesh bamg.geo bamg.mesh', shell=True)

    # delete bamg files for solution
    subprocess.call('rm *.bb', shell=True)
    # delete angener files
    subprocess.call('rm angener_prev.mesh profiles angener.mesh lines_file int_lines AMA_opers max_a_ratio.tri mesh', shell=True)
    # delete mmg files
    subprocess.call('rm *.sol adapted.mesh adapted.o.mesh output.msh output.mesh adj_hess_metric.mtr adj_metric.solb bamg.o.mesh', shell=True)
    # delete netgen files
    subprocess.call('rm netgen_prev.vol ng.ini ng.prof ngmesh.ini netgen.vol netgen.in2d netgen.geo', shell=True)  
    # delete metric related files
    subprocess.call('rm adj_metric.mtr implied_metric.mtr', shell=True)
    # delete solution and log files
    # subprocess.call('rm error-* log-*.txt', shell=True)
    # delete data dumps
    subprocess.call('rm error_data.csv int_err.csv cell_wise_error_*.csv max_ratios.csv solution_trial.txt mesh_data.csv error_data_solver.csv coefficients-*.txt', shell=True)
    # delete compilation files
    # subprocess.call('rm main.o main.so main.so.dSYM libpetsc_wrapper.a main.dylib ngs.prof.0', shell=True)
    # delete script files
    subprocess.call('cd Scripts; rm netgen2bamg netgen2angener netgen2mmg', shell=True)
    # delete some files from the cluster dump
    subprocess.call('rm core.cluster.rz.RWTH-Aachen.DE.* output_*.txt', shell=True)
    # delete some files from refine
    subprocess.call('rm adapted_prev.mesh adapted_prev.meshb adapted_surf.tec adapted.b8.ugrid adapted.meshb adapted_geom.tec adapted-final-metric.solb', shell=True)
    # delete some files for hp
    subprocess.call('rm nodeorder*.bb nodeorder*.sol nodeorder*.solb', shell=True)
    # delete existing solution files (if any)
    subprocess.call('rm solution-*', shell=True)
    # delete error csv files
    subprocess.call('rm cell_wise_*', shell=True)
    subprocess.call('rm ml_err_*', shell=True)
    subprocess.call('rm ngs.ini', shell=True)
    # delete log file
    subprocess.call('rm mesh-*', shell=True)
    subprocess.call('rm mldata-*', shell=True)

def replace(filename, string_to_search, replacement_string):
    with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace(string_to_search, replacement_string), end='')

def Calculate(no_steps):
    for p in range(0, no_steps):
        # subprocess.call('cd Scripts; g++ netgen2bamg_gmsh.cpp -o netgen2bamg_gmsh; ./netgen2bamg_gmsh', shell=True)
        # # call gmsh here
        # subprocess.call('python3 blossom.py', shell=True)

        # subprocess.call('cd Scripts; g++ bamg2netgen_quad.cpp -o bamg2netgen_quad; ./bamg2netgen_quad', shell=True)

        subprocess.call ("rm cell_wise_ml_err.csv", shell=True)
        
        
        
        # Testing for verification
        # subprocess.call ("ngs test.pde", shell=True)
        # subprocess.call('mv mldata-* ./ML/TrainingDataFiles', shell=True)
        
        
        # KUNAL 
        # Uncomment all next block to turn on the adaption
        subprocess.call ("cp netgen.in2d ./Scripts/netgen.in2d", shell=True)
        subprocess.call ("cp netgen.vol ./Scripts/netgen.vol", shell=True)
        subprocess.call ("cd ./Scripts; ./netgen2bamg_ho.out", shell=True)
        subprocess.call ("cp ./Scripts/bamg.mesh bamg.mesh", shell=True)
        subprocess.call ("cp ./Scripts/bamg.geo bamg.geo", shell=True)
        subprocess.call ("ngs test.pde", shell=True)
        if p < (no_steps - 1):
            subprocess.call (f"rm mldata-*", shell=True)
        else:
            subprocess.call('mv mldata-* ./ML/TrainingDataFiles', shell=True)
        
        
        
        
        
        
        
        
        
        
        
        
        
        # MakeNetgenMesh()
        
        
        
def Run(mesh_size, step, final, hp):
    no_calc = 1
    reset()
    if hp==0:
        for q in range(0, final):
            Calculate(no_calc)
            num_old = np.floor(mesh_size*(step)**q)
            num_new = np.floor(mesh_size*(step)**(q+1))
            x="define constant dof_target = {0}" .format(num_old)
            y="define constant dof_target = {0}" .format(num_new)

            replace("test.pde", x, y)
    if hp==1:
        Calculate(no_calc)
        x="define constant read_order = {0}" .format(0)
        y="define constant read_order = {0}" .format(1)
        replace("test.pde", x, y)
        for q in range(0, final):
            Calculate(no_calc)
            num_old = math.floor(mesh_size*(step)**q)
            num_new = math.floor(mesh_size*(step)**(q+1))
            x="define constant dof_target = {0}" .format(num_old)
            y="define constant dof_target = {0}" .format(num_new)
            replace("test.pde", x, y)

    x="define constant read_order = {0}" .format(1)
    y="define constant read_order = {0}" .format(0)
    replace("test.pde", x, y)

    num_final = math.floor(mesh_size*(step)**(final))
    num_initial = math.floor(mesh_size)
    x="define constant dof_target = {0}" .format(num_final)
    y="define constant dof_target = {0}" .format(num_initial)
    replace("test.pde", x, y)


def reset():
    DeleteFiles()
    # subprocess.call('cp Meshes/netgen/square32.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/square.in2d netgen.in2d', shell=True)
    # subprocess.call('cp Meshes/netgen/square_large32.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/square_large.in2d netgen.in2d', shell=True)
    # subprocess.call('cp Meshes/netgen/flatplate.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/flatplate.in2d netgen.in2d', shell=True)
    subprocess.call('cp Meshes/netgen/naca_georg_1000.vol netgen.vol', shell=True)
    subprocess.call('cp Meshes/netgen/naca_georg_1000.in2d netgen.in2d', shell=True)

    # subprocess.call('cp Meshes/netgen/diamond.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/diamond.in2d netgen.in2d', shell=True)
    # subprocess.call('cp Meshes/netgen/cube96.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/cube.geo netgen.geo', shell=True)
    # subprocess.call('cp Meshes/netgen/delta.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/refine/delta01.meshb adapted.meshb', shell=True)
    # subprocess.call("cp Meshes/refine/delta01.meshb adapted.meshb", shell=True)
    # subprocess.call("transmesh adapted.meshb adapted.mesh", shell=True)
    # subprocess.call('cd Scripts; g++ refine2netgen.cpp -o refine2netgen; ./refine2netgen', shell=True)
    # subprocess.call("cp netgen.vol delta_wing.vol", shell=True)
    # subprocess.call('python3 delta_fix.py', shell=True)
    # subprocess.call('mv delta_wing_fix.vol netgen.vol', shell=True)
    # subprocess.call('rm delta_wing.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/wedge.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/wedge.in2d netgen.in2d', shell=True)

    subprocess.call('cd Scripts; g++ netgen2bamg.cpp -o netgen2bamg; ./netgen2bamg', shell=True)
    # subprocess.call('cd Scripts; g++ netgen2mmg2d.cpp -o netgen2mmg2d; ./netgen2mmg2d', shell=True)
    # subprocess.call('cd Scripts; g++ netgen2mmg3d.cpp -o netgen2mmg3d; ./netgen2mmg3d', shell=True)
    # subprocess.call('cd Scripts; g++ netgen2madlib2d.cpp -o netgen2madlib2d; ./netgen2madlib2d', shell=True)
    # subprocess.call('cp Meshes/netgen/naca_georg_coarse.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/naca_georg_coarse.in2d netgen.in2d', shell=True)
    # subprocess.call('cp Meshes/netgen/bump.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/bump.in2d netgen.in2d', shell=True)
    # subprocess.call('cp Meshes/netgen/naca_georg_1000.vol netgen.vol', shell=True)
    # subprocess.call('cp Meshes/netgen/naca_georg_1000.in2d netgen.in2d', shell=True)

    # # # subprocess.call('cp Meshes/netgen/bump_ho.vol netgen.vol', shell=True)
    # # # subprocess.call('cp Meshes/netgen/bump_ho.in2d netgen.in2d', shell=True)
    # # # # subprocess.call('cp Meshes/netgen/cyl_full.vol netgen.vol', shell=True)
    # # # # subprocess.call('cp Meshes/netgen/cyl_full_quad.in2d netgen.in2d', shell=True)
    # subprocess.call('cd Scripts; g++ netgen2bamg_ho.cpp -o netgen2bamg_ho; ./netgen2bamg_ho', shell=True)
    # subprocess.call('bamg -g bamg.geo -o bamg.mesh -coef 4.28', shell=True)
    # subprocess.call('cd Scripts; g++ bamg2netgen.cpp -o bamg2netgen; ./bamg2netgen', shell=True)
    # subprocess.call('cd Scripts; g++ netgen2madlib3d.cpp -o netgen2madlib3d; ./netgen2madlib3d', shell=True)
    subprocess.call('make clean; make', shell=True)

def StandardMeshes():
    DeleteFiles()
    # subprocess.call('cp Meshes/netgen/square.in2d netgen.in2d', shell=True)
    subprocess.call('cp Meshes/netgen/naca_georg_1000.in2d netgen.in2d', shell=True)
    # subprocess.call('cp Meshes/netgen/flatplate.in2d netgen.in2d', shell=True)
    subprocess.call('make clean; make', shell=True)
    mesh_sizes = ["2", "8", "32", "128", "512", "2048"]
    # mesh_sizes = ["35x25", "69x49", "137x97", "273x193", "545x385"]
    for q in range(0, len(mesh_sizes)):
        subprocess.call('cp Meshes/netgen/naca_georg_1000'+mesh_sizes[q]+'.vol netgen.vol', shell=True)
        Calculate(1)

def GenerateMLData():
    DeleteFiles()

    subprocess.call('make clean; make', shell=True)

    subprocess.call('cp Meshes/netgen/visc_subsonic_naca.in2d netgen.in2d', shell=True)

    mesh_sizes = ["2334"]
    
    N = 10
    
    mach_range = np.linspace(0.35, 0.7, N)
    aoa_range = np.linspace(0.0, 6.0, N)

    for q in range(0, len(mesh_sizes)):
        # Write the loops over parameters
        for i in range(0, len(mach_range)):
            for j in range(0, len(aoa_range)):
                subprocess.call('cp Meshes/netgen/visc_subsonic_naca_'+mesh_sizes[q]+'.vol netgen.vol', shell=True)
                subprocess.call('cp ML/test_ml.pde test.pde', shell=True)
                mach_old = 0.35
                mach_new = mach_range[i]

                x="define constant mach = {0}" .format(mach_old)
                y="define constant mach = {0}" .format(mach_new)
                replace("test.pde", x, y)

                aoa_old = 0.0
                aoa_new = aoa_range[j]

                x="define constant alpha = {0}" .format(aoa_old)
                y="define constant alpha = {0}" .format(aoa_new)
                replace("test.pde", x, y)

                Calculate(5)
                subprocess.call('rm mesh*', shell=True)
                # subprocess.call('rm log*', shell=True)
                # KUNAL
                # We do not need mesh file to train the model.
                # Uncomment the below line to see the end mesh
                # subprocess.call(f'cp bamg.mesh bamg_eps_{eps_range[i]}_beta1_{beta1_range[j]}_beta2_{beta2_range[k]}.mesh', shell=True)

    # KUNAL
    # subprocess.call('rm ./ML/TrainingDataFiles/*', shell=True)
    # subprocess.call('rmdir ./ML/TrainingDataFiles', shell=True)
    # subprocess.call('mkdir ./ML/TrainingDataFiles', shell=True)
    # subprocess.call('mv mldata-* ./ML/TrainingDataFiles', shell=True)
    
    # KUNAL
    # subprocess.call('rm ./ML/LogDataFiles/*', shell=True)
    # subprocess.call('rmdir ./ML/LogDataFiles', shell=True)
    # subprocess.call('mkdir ./ML/LogDataFiles', shell=True)
    # subprocess.call('mv log-* ./ML/LogDataFiles', shell=True)
    
    print(f"Neglect if there is some rm or rmdir error.")

# KUNAL
    
def Processing_Data():
    print("Removing any previous input files: *.npy")
    subprocess.call('cd ML; rm *.npy', shell=True)
    print(f"Neglect if there is some rm or rmdir error.")
    
    print("Creating the input and output files for the ML model: Initiated")
    subprocess.call('cd ML; python3 file_reader.py', shell=True)
    print("Creating the input and output files for the ML model: Completed")
    

def TrainMLModel():
    
    print("Removing any previous weight files: *.h5")
    subprocess.call('rm *.h5', shell=True)
    print(f"Neglect if there is some rm or rmdir error.")
    
    print("Training the ML model: Initiated")
    subprocess.call('python3 training.py', shell=True)
    subprocess.call('mv *.h5 ./ML/Multiple_Models_Weights/',shell=True)
    subprocess.call('mv ./ML/*.save ./ML/Multiple_Models_Output_Scaler/',shell=True)
    print("Training the ML model: Completed")    

    
def DeployML():
    DeleteFiles()

    print(f"Removing all the old files of in this directory\n")
    subprocess.call('rm cell_wise_ml_err_*',shell=True)
    subprocess.call('rm error_data.csv*',shell=True)
    subprocess.call('rm ml_err_data.csv*',shell=True)

    subprocess.call('make clean; make', shell=True)

    subprocess.call('cp Meshes/netgen/square.in2d netgen.in2d', shell=True)

    
    final_size = [128]
    # final_size = [32]
    for ne in final_size:
    
    	# KUNAL 
        # I am commenting the next line
        subprocess.call('cp Meshes/netgen/square'+str(ne)+'.vol netgen.vol', shell=True)
        
        # KUNAL
        subprocess.call('cd Scripts; g++ netgen2bamg_ho.cpp -o netgen2bamg_ho', shell=True)
        
        # Write the loops over parameters
        Calculate(10)
        subprocess.call("cp cell_wise_ml_err.csv cell_wise_ml_err_"+str(ne)+".csv", shell=True)
        subprocess.call("rm cell_wise_ml_err.csv", shell=True)
        
        
        # subprocess.call('cp adj_metric.mtr adj_metric_'+str(ne)+'.mtr', shell=True)
        
        # subprocess.call('cp netgen.vol ./Scripts/', shell=True)
        # subprocess.call('cd Scripts; ./netgen2bamg', shell=True)
        # subprocess.call('cd Scripts; cp bamg.mesh  ../', shell=True)
        # subprocess.call('cd Scripts; cp bamg.geo  ../', shell=True)
        subprocess.call('cp bamg_prev.mesh bamg_prev_'+str(ne)+'.mesh', shell=True)
        
def post_processing():

    print("Creating the NDOF VS Error plot")
    
    subprocess.call('cd Post_Processing_Scripts; cp cell_wise_ml_error_all.py ../',shell=True)

    subprocess.call('python3 cell_wise_ml_error_all.py',shell=True)

    subprocess.call('rm cell_wise_ml_error_all.py',shell=True)
    
    print("Creating the Element wise ml error plot")
    
    subprocess.call('cd Post_Processing_Scripts; cp NDOF_Error_Plot.py ../',shell=True)

    subprocess.call('python3 NDOF_Error_Plot.py',shell=True)
    
    subprocess.call('rm NDOF_Error_Plot.py',shell=True)

    # subprocess.call('cp *.png ../',shell=True)

def post_processing_clean():
    print("Removing all the graphs in this folder.\n")
    subprocess.call('rm *.png',shell=True)
    


# DeleteFiles()
GenerateMLData()
# KUNAL
# Processing_Data()
# TrainMLModel()
# subprocess.call('python3 Deploy_File_Changes.py',shell=True)
# DeployML()
# post_processing()
# post_processing_clean()

end_time = time.time()
elapsed_time = end_time-start_time
print(f"Elapsed time: {elapsed_time} seconds")
