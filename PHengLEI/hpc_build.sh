#!/bin/bash
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
#          P   P  H   H  E      NN   N  G      L      E       I          +
#          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
#          P      H   H  E      N  N N  G   G  L      E       I          +
#          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
#------------------------------------------------------------------------+
#          Platform for Hybrid Engineering Simulation of Flows           +
#           China Aerodynamics Research and Development Center           +
#                     (C) Copyright, Since 2010                          +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# File: hpc_build.sh
# Purpose: Supercomputing Environment Deployment and Build Automation
#! @author    Luo Junyi.

ENV_FILE_NAME="env.sh"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ENV_FILE="$SCRIPT_DIR/$ENV_FILE_NAME"

# Platform Selection Interface
echo "[PHengLEI Supercomputing Platform Deployment Process]"
echo "Select HPC Platform Configuration:"
echo "1. NSCC-CD, X86, GCC    (Chengdu Supercomputing Center, x86_64, GNU Compiler Collection)"
echo "2. NSCC-CD, X86, INTEL  (Chengdu Supercomputing Center, x86_64, Intel Classic Compilers)"
echo "3. NSCC-TJ, ARM, GCC    (Tianjin Supercomputing Center, ARM64, GNU Compiler Collection)"
echo "4. Tianhe-2, X86, GCC   (Guangzhou Tianhe-2, x86_64, GNU Compiler Collection)"
echo "5. NSCC-WX, X86, GCC    (Wuxi Supercomputing Center, x86_64, GNU Compiler Collection)"
echo "6. NSCC-JN, X86, INTEL  (Jinan Supercomputing Center, x86_64, Intel Classic Compilers)"

# User Input Handling
read -p "Enter platform selection (1-6): " platform

# Initialize Environment Configuration File
echo "#!/bin/bash" > $ENV_FILE
echo "# Auto-generated HPC Environment Configuration" >> $ENV_FILE

# Platform-Specific Module Configuration
case $platform in
    1)
        echo "Configuring NSCC-CD GCC environment..."
        {
            echo "module purge"
            echo "module load compiler/devtoolset/7.3.1"
            echo "module load mpi/hpcx/2.4.1/gcc-7.3.1"
            echo "module load compiler/cmake/3.20.1"
        } >> $ENV_FILE
        ;;
    2)
        echo "Configuring NSCC-CD Intel environment..."
        {
            echo "module purge"
            echo "module load compiler/intel/2017.5.239"
            echo "module load mpi/hpcx/2.11.0/intel-2017.5.239"
            echo "module load compiler/cmake/3.20.1"
        } >> $ENV_FILE
        ;;
    3)
        echo "Configuring NSCC-TJ ARM environment..."
        {
            echo "module purge"
            echo "module load GCC/9.3.0"
            echo "module load mpich/mpi-x-gcc9.3.0"
            echo "module load cmake/3.20.4"
        } >> $ENV_FILE
        ;;
    4)
        echo "Configuring Tianhe-2 GCC environment..."
        {
            echo "module purge"
            echo "module load gcc/9.5.0"
            echo "module load mpi/mpich/4.1.2-gcc-9.5.0-ch3"
            echo "module load cmake/3.30.1-gcc-9.5.0"
        } >> $ENV_FILE
        ;;
    5)
        echo "Configuring Wuxi Supercomputer environment..."
        {
            echo "spack unload --all"
            echo "spack load mpich@3.4.2%gcc@=4.8.5 arch=linux-centos7-haswell"
            echo "spack load cmake@3.22.1%gcc@4.8.5 arch=linux-centos7-haswell"
        } >> $ENV_FILE
        ;;
    6)
        echo "Configuring NSCC-JN Intel environment..."
        {
            echo "module purge"
            echo "module load intel/2018"
            echo "module load intelmpi/2018"
            echo "module load cmake/3.21.1"
        } >> $ENV_FILE
        ;;
    *)
        echo "Error: Invalid platform selection"
        rm -f $ENV_FILE
        exit 1
        ;;
esac

# Environment Verification Utilities
{
    echo "echo \"Current Computational Environment Status:\""
    echo "echo -n \"• GCC Version:    \"; gcc --version | head -n1"
    echo "echo -n \"• MPI Version:    \"; mpicc -v 2>&1 | head -n1"
    echo "echo -n \"• CMake Version:  \"; cmake --version | head -n1"
} >> $ENV_FILE

# Permission Management and User Guidance
chmod +x $ENV_FILE
echo "Environment configuration archived to: $ENV_FILE_NAME"
echo "To restore this compilation environment subsequently, execute:"
echo "  $ source $ENV_FILE_NAME"

# User confirmation
echo "[Press Enter to continue or wait 15 seconds]"
read -t 15 dummy || true

# Environment Activation and Build Execution
echo -e "\nInitializing modules and commencing build process..."
. "$ENV_FILE" && make -C "$SCRIPT_DIR" config hpc_mode=ON
