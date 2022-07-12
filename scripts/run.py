import os
import time
import numpy as np
import os
import math
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties

# basic configuration 
project_path = "/rds/general/user/dx121/home/mpi-coursework-acse-dx121"
# project_path = "/Users/xingzheng/mpi-coursework-acse-dx121"
max_processors = 62
max_size = 550

# collect_running_time run a user-defined command and return the run time.
def Collect_Runing_time(run_command):
    start = time.time()
    res = os.popen(run_command).read()
    print(res)
    
    return time.time() - start

# plot is used to plot the run time.
def plot(times_array, size):
    fig = plt.figure(figsize=(20, 15))
    ax1 = fig.add_subplot(111)
    ax1.margins(0.1)
    ax1.grid(True)
    
    ax1.plot([i for i in range(1, max_processors)], times_array, 'k.-', label='parallel wave equation - basic vesion')
    
    ax1.set_xlabel('the number of processor', fontsize=16)
    ax1.set_ylabel('$time(second)$', fontsize=16)
    
    ax1.set_title('time cost of parallel wave equation', fontsize=16)
    ax1.legend(loc='best', fontsize=14)
    plt.savefig(project_path + "/img/parallel_wave_equation_basic_version.png")
    
def create_script(processor, size):
    scrpts = "#PBS -N 2-mpi_%d_%d\n#PBS -l walltime=1:00:00\n#PBS -l select=2:ncpus=%d:mpiprocs=%d:mem=%dGB\nmodule load intel-suite/2019.4\nmodule load mpi/intel-2019.8.254\nmodule load cmake/3.18.2 \ncd /rds/general/user/dx121/home/mpi-coursework-acse-dx121\ncd build\ncmake .. && make && cp MPI-CourseWork ../hpc-work/\ncd ../hpc-work/\nmpiexec ~/mpi-coursework-acse-dx121/hpc-work/MPI-CourseWork %d"
    filename = "../hpc-work/2-run-%d-%d.pbs" %(size, processor)
    with open(filename, "w") as text_file:
        text_file.write(scrpts % (size, processor, processor, processor, 8, size))
    return filename


def main():
    # for size in range(200, max_size, 50):
    #     times = []
    #     for p in range(1, max_processors):
    #         command = 'mpiexec -n %d %s/hpc-work/MPI-CourseWork %d'%(p, project_path, size)
    #         times.append(Collect_Runing_time(command))
        
    #     times_array = np.array(times)
    #     print("Time collected: ", times_array, "   size:", size, "   max_processors:", max_processors)
    #     np.save(project_path + "/data/basic-%d-%d.npy"%(max_processors, size), times_array)
    #     # plot(times_array, size)
    os.chdir("../hpc-work")
    size = [100, 200, 300, 400, 500]
    # processor = [1, 2, 4, 6, 8, 16, 32]
    processor = [32]
    for s in size:
        for p in processor:
            filename = create_script(p, s)
            print(filename)
            Collect_Runing_time("qsub %s"%filename)
            # os.remove(filename)
        
    
if __name__ == "__main__" :
    main()
    