# Only focus on the first three lines

# To select the compilation option
# "None"  -> Classical compilation
# "Debug" -> Debug mode compilation
# "Optim" -> Optimization compilation
version="Debug"

# To enable or not the use of the MUMPS Library
# 1 -> With
# 0 -> Without
USE_MUMPS=1

# If MUMPS is employed, provide path to the home directory
MUMPS_DIR="/home/tiffanie/library/MUMPS_4.10.0_SEQ"
BLIBS = "L/usr/local/lib -lblas -llapack"


####################################################
####################################################




if [ $version != "None" ] && [ $version != "Optim" ] && [ $version != "Debug" ] ; then
    echo " "
    echo " ----------------------------------------------------"
    echo " ERROR in the selected option for compilation mode..."
    echo "   Please modify and choose between"
    echo '     "None"'
    echo '     "Optim"'
    echo '     "Debug"'
    echo " "
    exit
fi

if [ $USE_MUMPS != 0 ] && [ $USE_MUMPS != 1 ] ; then
    echo " "
    echo " ----------------------------------------------------"
    echo " ERROR in the selected option for MUMPS..."
    echo "   Please modify and choose between"
    echo '     1 -> Compilation with MUMPS'
    echo '     2 -> Compilation without MUMPS'
    echo " ----------------------------------------------------"
    echo " "
    exit
fi  

echo " "
echo " ------------------------------------------------------------------------------------------"
echo "generation of the Makefile for DaFFEM using CMake"
echo " "
echo " -- Chosen optimisation option:" $version
if [ $USE_MUMPS =  1 ] ; then
    echo " -- Compilation is done with MUMPS."
    echo " -- MUMPS home path given :" $MUMPS_DIR
else
    echo " -- Compilation without MUMPS."
    echo "Be careful to adjust datafile with GMRES-ILU before running"
fi
echo " "
echo " ------------------------------------------------------------------------------------------"
echo " "


if [ $version = "None" ] ; then
    flag=" "
elif [ $version = "Optim" ] ; then
    flag="-O3"
elif [ $version = "Debug" ] ; then
    flag="-g -fbacktrace -ffpe-trap=zero,overflow,underflow" 
fi

if [ $USE_MUMPS = 1 ] ; then
    cmake .. -DUSE_MUMPS=$USE_MUMPS -DMUMPS_DIR=$MUMPS_DIR -Dadd_flag=$add_flag
else
    cmake .. -DUSE_MUMPS=$USE_MUMPS -Dadd_flag=$add_flag  
fi
