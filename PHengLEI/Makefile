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
#! @file      Makefile
#! @brief     It defines the configuration options.
#! @author    Bell, Luo Junyi.


# Compilation Configuration Options for the PHengLEI Solver
BUILD_API ?= ON
BUILD_CFD ?= ON
BUILD_DemoSolver ?= OFF
BUILD_GMRESSolver ?= OFF
BUILD_GPU ?= OFF
BUILD_IncomSolver ?= ON
BUILD_LBMSolverMPI ?= OFF
BUILD_LBMSolverOMP ?= OFF
BUILD_LIBRARY_TYPE ?= Static
BUILD_MESH ?= ON
BUILD_SpecDiffHybSolver ?= OFF
CMAKE_BUILD_TYPE ?= Release
USE_AI_RandomForestRegressor ?= OFF
USE_DOUBLE_PRECISION ?= ON
USE_ERROR_DATA_DEBUG ?= OFF
USE_HypreLib ?= ON
USE_OCC ?= OFF
USE_OMP ?= OFF
USE_SMARTARRAY ?= ON
USE_TecplotLib ?= ON
USE_UnapLib ?= OFF
USE_Windows_X64 ?= OFF
USE_YhamgLib ?= OFF

PHengLEI_CMAKE_OPTIONS = \
	-DBUILD_API=$(BUILD_API) \
	-DBUILD_CFD=$(BUILD_CFD) \
	-DBUILD_DemoSolver=$(BUILD_DemoSolver) \
	-DBUILD_GMRESSolver=$(BUILD_GMRESSolver) \
	-DBUILD_GPU=$(BUILD_GPU) \
	-DBUILD_IncomSolver=$(BUILD_IncomSolver) \
	-DBUILD_LBMSolverMPI=$(BUILD_LBMSolverMPI) \
	-DBUILD_LBMSolverOMP=$(BUILD_LBMSolverOMP) \
	-DBUILD_LIBRARY_TYPE=$(BUILD_LIBRARY_TYPE) \
	-DBUILD_MESH=$(BUILD_MESH) \
	-DBUILD_SpecDiffHybSolver=$(BUILD_SpecDiffHybSolver) \
	-DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) \
	-DUSE_AI_RandomForestRegressor=$(USE_AI_RandomForestRegressor) \
	-DUSE_DOUBLE_PRECISION=$(USE_DOUBLE_PRECISION) \
	-DUSE_ERROR_DATA_DEBUG=$(USE_ERROR_DATA_DEBUG) \
	-DUSE_HypreLib=$(USE_HypreLib) \
	-DUSE_OCC=$(USE_OCC) \
	-DUSE_OMP=$(USE_OMP) \
	-DUSE_SMARTARRAY=$(USE_SMARTARRAY) \
	-DUSE_TecplotLib=$(USE_TecplotLib) \
	-DUSE_UnapLib=$(USE_UnapLib) \
	-DUSE_Windows_X64=$(USE_Windows_X64) \
	-DUSE_YhamgLib=$(USE_YhamgLib)


# 3rdparty lib version info
hdf5_version     = 1.14.4-3
cgns_version     = 4.4.0
parmetis_version = 4.0.3
metis_version    = 5.0.1
mgird_version    = 1.0
hypre_version    = 2.24.0

# Ensure bash shell for advanced features
SHELL := /bin/bash

# Number of processor to parallel compiling.
nProcessor   = $(shell nproc)
cc           = mpicc
cxx          = mpicxx

# Environment checks configurations
REQUIRED_GCC_VERSION := 4.9.2
RECOMMENDED_MPI_TYPE  := MPICH
RECOMMENDED_MPI_VERSION := 3.3.2

# Basically envioment informations.
cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)
PATH_MPICC           = $(shell which mpicc)
PATH_MPICXX          = $(shell which mpicxx)
PATH_MPIEXEC         = $(shell which mpiexec)

# Path for 3rdparty.
PARMETISBUILDDIR     = build/$(systype)-$(cputype)
LIB_SOURCE_DIR       = 3rdparty
LIB_BINARY_DIR       = $(LIB_SOURCE_DIR)/liblinux/$(systype)-$(cputype)
LIB_HDF5_INCLUDE_DIR = $(LIB_SOURCE_DIR)/HDF5/linux/$(systype)-$(cputype)/gcc/
LIB_CGNS_INCLUDE_DIR = $(LIB_SOURCE_DIR)/cgns/linux/$(systype)-$(cputype)/gcc/

# Architecture detection and info output.
ifeq ($(cputype),aarch64)
LIB_BINARY_DIR       = $(LIB_SOURCE_DIR)/liblinux/$(systype)-$(cputype)
$(info -------------------------------------------------------------------------------------)
$(info [Build Environment Configuration Info]                                               )
$(info • System Type:          $(systype)                                                   )
$(info • Processor Architecture: ARM64 (aarch64)                                            )
$(info • Auto-configured 3rdparty Library Path: $(LIB_BINARY_DIR)                           )
$(info -------------------------------------------------------------------------------------)
else ifeq ($(cputype),x86_64)
# Compiler detection, whether GNU Compiler Collection or Intel Classic Compilers.
MPI_PATH_LOWER := $(shell echo "$(PATH_MPICC)$(PATH_MPICXX)" | tr '[:upper:]' '[:lower:]')
ifneq (,$(findstring intel,$(MPI_PATH_LOWER)))
COMPILER_TYPE := Intel
LIB_BINARY_DIR = $(LIB_SOURCE_DIR)/liblinux/$(systype)-$(cputype)/intel/
LIB_HDF5_INCLUDE_DIR = $(LIB_SOURCE_DIR)/HDF5/linux/$(systype)-$(cputype)/intel/
LIB_CGNS_INCLUDE_DIR = $(LIB_SOURCE_DIR)/cgns/linux/$(systype)-$(cputype)/intel/
CMAKE_INTEL_FLAG = -DUSE_INTEL_COMPILER=ON
else
COMPILER_TYPE := GCC
LIB_BINARY_DIR = $(LIB_SOURCE_DIR)/liblinux/$(systype)-$(cputype)/gcc/
CMAKE_INTEL_FLAG =
endif
$(info -------------------------------------------------------------------------------------)
$(info [Build Environment Configuration Info]                                               )
$(info • System Type:          $(systype)                                                   )
$(info • Processor Architecture: x86_64                                                     )
$(info • MPI Compiler Type:    $(COMPILER_TYPE)                                             )
$(info • MPI Compiler Path:    $(PATH_MPICC)                                                )
$(info • Auto-configured Library Path: $(LIB_BINARY_DIR)                                    )
$(info -------------------------------------------------------------------------------------)

else
$(error Unsupported CPU architecture: $(cputype))
endif

# --------------------- Environment Validation ---------------------
define check-env
	@# Required tools existence check
	@echo "Checking build environment prerequisites..."
	@command -v gcc >/dev/null 2>&1 || { echo "Error: Required compiler 'gcc' not found. Build aborted."; exit 1; }
	@command -v cmake >/dev/null 2>&1 || { echo "Error: Required build tool 'cmake' not found. Build aborted."; exit 1; }
	@command -v mpicc >/dev/null 2>&1 || { echo "Error: MPI compiler 'mpicc' not found. Build aborted."; exit 1; }
	@command -v mpicxx >/dev/null 2>&1 || { echo "Error: MPI compiler 'mpicxx' not found. Build aborted."; exit 1; }

	@# CMake version check
	@echo "Verifying CMake version..."
	@current_cmake=$$(cmake --version | head -n1 | awk '{print $$3}'); \
	required_cmake="3.18"; \
	if [ "$$(printf '%s\n' "$$required_cmake" "$$current_cmake" | sort -V | head -n1)" != "$$required_cmake" ]; then \
		echo "Error: Detected CMake version $${current_cmake}"; \
		echo "Minimum required CMake version: $${required_cmake}"; \
		read -p "Press any key to acknowledge this error and exit..." dummy; \
		exit 1; \
	fi

	@# GCC version validation
	@echo "Verifying GCC version..."
	@current_gcc=$$(gcc -dumpversion); \
	if [ "$$(printf '%s\n' "$(REQUIRED_GCC_VERSION)" "$$current_gcc" | sort -V | head -n1)" != "$(REQUIRED_GCC_VERSION)" ]; then \
		echo "Warning: Detected GCC version $${current_gcc} (minimum required: $(REQUIRED_GCC_VERSION))"; \
		read -p "Potential compilation risks exist. Continue? [Y/N]: " -t 30 confirm; \
		case "$$confirm" in \
			[Yy]*) ;; \
			[Nn]*) echo "Build aborted by user."; exit 1;; \
			*) echo "Invalid input, aborting build."; exit 1;; \
		esac; \
    fi

    @# MPI implementation check
    @echo "Verifying mpi version..."
    @mpi_version_output=$$(mpicc -v 2>&1); \
    mpi_vendor=$$(echo "$$mpi_version_output" | grep -oi 'MPICH\|Open MPI' | head -n1); \
    if [ -z "$$mpi_vendor" ]; then \
        mpi_vendor="Unrecognized"; \
    fi; \
    mpi_version=$$(echo "$$mpi_version_output" | grep -i 'MPICH version' | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1); \
    \
    if [ "$$mpi_vendor" != "$(RECOMMENDED_MPI_TYPE)" ]; then \
        if [ "$$mpi_vendor" = "Unrecognized" ]; then \
            echo "Warning: Unrecognized MPI implementation detected"; \
            echo "Required MPI: $(RECOMMENDED_MPI_TYPE)-$(RECOMMENDED_MPI_VERSION)"; \
            read -p "Continue with UNIDENTIFIED MPI? [Y/N]: " -t 30 confirm; \
        else \
            echo "Warning: Detected $$mpi_vendor (Requires $(RECOMMENDED_MPI_TYPE))"; \
            read -p "Unsupported MPI implementation. Continue? [Y/N]: " -t 30 confirm; \
        fi; \
        case "$$confirm" in \
            [Yy]*) echo "Proceeding with $$mpi_vendor implementation";; \
            *) echo "Build aborted."; exit 1;; \
        esac; \
    elif [ -z "$$mpi_version" ]; then \
        echo "Warning: Failed to detect $(RECOMMENDED_MPI_TYPE) version"; \
        read -p "Version verification unavailable. Continue? [Y/N]: " -t 30 confirm; \
        case "$$confirm" in \
            [Yy]*) echo "Proceeding with UNVERIFIED MPI implementation";; \
            *) echo "Build aborted by user."; exit 1;; \
        esac; \
    elif [ "$$mpi_version" != "$(RECOMMENDED_MPI_VERSION)" ]; then \
        echo "Warning: $(RECOMMENDED_MPI_TYPE) version mismatch (Current: $${mpi_version}, Required: $(RECOMMENDED_MPI_VERSION))"; \
        read -p "Continue with version $${mpi_version}? [Y/N]: " -t 30 confirm; \
        case "$$confirm" in \
            [Yy]*) ;; \
            *) echo "Build aborted."; exit 1;; \
        esac; \
    fi

endef

# --------------------- Build Configuration ------------------------
define run-config
	@# Important: the default path is the one of Makefile's, except for '&&' connectivity.
	mkdir -p $(LIB_BINARY_DIR)
	mkdir -p $(LIB_HDF5_INCLUDE_DIR)
	mkdir -p $(LIB_CGNS_INCLUDE_DIR)

	@# Boost/blitz/eigen
	@echo "Extract eigen library..."
	cd $(LIB_SOURCE_DIR) && tar -xvf eigen.tar.gz

	@# Tecplot.
	@# Tecplot lib could not been compiled, must be installed for different linux platform, so copy it directly.
	@# it's supported only for Linux-aarch64 and Linux-x86_64 now.
	@# cd $(LIB_SOURCE_DIR) && cp libtecio-$(systype)-$(cputype).a $(CURDIR)/$(LIB_BINARY_DIR)/libtecio.a

	@# HDF5
	@echo "Build HDF5 library..."
	cd $(LIB_SOURCE_DIR) && unzip -o hdf5-$(hdf5_version).zip
	cd $(LIB_SOURCE_DIR)/hdf5-$(hdf5_version) && mkdir build && cd build && mkdir install && cmake .. -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_C_FLAGS="-ldl" -DBUILD_SHARED_LIBS=OFF -DHDF5_BUILD_HL_LIB=ON -DHDF5_BUILD_STATIC_TOOLS=ON -DHDF5_BUILD_TOOLS=OFF -DHDF5_ENABLE_PARALLEL=ON -DHDF5_ENABLE_SZIP_ENCODING=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=OFF -DHDF_BUILD_HL_LIB=ON -DHDF_ENABLE_PARALLEL=ON -DUSE_SHARED_LIBS=OFF && make -j $(nProcessor) && make install  
	@cp $(CURDIR)/$(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)/build/install/lib/*.a $(CURDIR)/$(LIB_BINARY_DIR)
	@cp -r $(CURDIR)/$(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)/build/install/include/ $(CURDIR)/$(LIB_HDF5_INCLUDE_DIR) 

	@# CGNS
	@echo "Build CGNS library..."
	cd $(LIB_SOURCE_DIR) && unzip -o CGNS-$(cgns_version).zip
	cd $(LIB_SOURCE_DIR)/CGNS-$(cgns_version) && mkdir build && cd build && mkdir install && cmake .. -DCMAKE_INSTALL_PREFIX=./install -DCGNS_ENABLE_HDF5=ON -DHDF5_DIR=$(CURDIR)/$(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)/build/install/cmake -DHDF5_NEED_MPI=ON -DCGNS_BUILD_SHARED=OFF -DHDF5_NEED_MPI=ON -DHDF5_IS_PARALLEL=ON && make -j $(nProcessor) && make install 
	@cp $(CURDIR)/$(LIB_SOURCE_DIR)/CGNS-$(cgns_version)/build/install/lib/*.a $(CURDIR)/$(LIB_BINARY_DIR)
	@cp -r $(CURDIR)/$(LIB_SOURCE_DIR)/CGNS-$(cgns_version)/build/install/include/ $(CURDIR)/$(LIB_CGNS_INCLUDE_DIR)

	@# Metis && Parmetis.
	@echo "Build Metis and Parmetis library..."
	cd $(LIB_SOURCE_DIR) && tar -xvf parmetis-$(parmetis_version).tar.gz
	cd $(LIB_SOURCE_DIR)/parmetis-$(parmetis_version) && make config prefix=$(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)
	cd $(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)/$(PARMETISBUILDDIR) && make -j $(nProcessor)  
	@cp $(CURDIR)/$(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)/$(PARMETISBUILDDIR)/libparmetis/*.a $(CURDIR)/$(LIB_BINARY_DIR)
	@cp -r $(CURDIR)/$(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)/include/ $(CURDIR)/$(LIB_SOURCE_DIR)/parmetis
	cd $(LIB_SOURCE_DIR) && tar -xvf metis-$(metis_version).tar.gz
	cd $(LIB_SOURCE_DIR)/metis-$(metis_version) && make config prefix=$(LIB_SOURCE_DIR)/metis-$(metis_version)
	cd $(LIB_SOURCE_DIR)/metis-$(metis_version)/$(PARMETISBUILDDIR) && make -j $(nProcessor)
	@cp $(CURDIR)/$(LIB_SOURCE_DIR)/metis-$(metis_version)/$(PARMETISBUILDDIR)/libmetis/*.a $(CURDIR)/$(LIB_BINARY_DIR)
	@cp -r $(CURDIR)/$(LIB_SOURCE_DIR)/metis-$(metis_version)/include/ $(CURDIR)/$(LIB_SOURCE_DIR)/metis

	@# mgrid
	@echo "Build MGridGen library..."
	cd $(LIB_SOURCE_DIR) && tar -xvf mgrid-$(mgird_version).tar.gz
	cd $(LIB_SOURCE_DIR)/mgrid-$(mgird_version)/build && cmake .. && make -j $(nProcessor)
	@cp $(CURDIR)/$(LIB_SOURCE_DIR)/mgrid-$(mgird_version)/build/MGridGen/Lib/*.a $(CURDIR)/$(LIB_BINARY_DIR)
	@cp $(CURDIR)/$(LIB_SOURCE_DIR)/mgrid-$(mgird_version)/build/MGridGen/IMlib/*.a $(CURDIR)/$(LIB_BINARY_DIR)

	@# hypre
	@echo "Build Hypre library..."
	cd $(LIB_SOURCE_DIR) && unzip -o hypre-$(hypre_version).zip
	cd $(LIB_SOURCE_DIR)/hypre-$(hypre_version)/src/cmbuild && cmake .. && make -j $(nProcessor) && make install
	@# detect hypre lib path dynamically
	@HYPRE_LIB_PATH=$$( \
	if [ -d "$(LIB_SOURCE_DIR)/hypre-$(hypre_version)/src/hypre/lib64" ]; then echo "$(LIB_SOURCE_DIR)/hypre-$(hypre_version)/src/hypre/lib64"; \
	elif [ -d "$(LIB_SOURCE_DIR)/hypre-$(hypre_version)/src/hypre/lib" ]; then echo "$(LIB_SOURCE_DIR)/hypre-$(hypre_version)/src/hypre/lib"; \
	else echo "ERROR: Hypre library path not found"; exit 1; fi); \
	echo "Detected hypre library path: $$HYPRE_LIB_PATH"; \
	cp $$HYPRE_LIB_PATH/libHYPRE.a $(CURDIR)/$(LIB_BINARY_DIR)/libhypre.a && \
	cp -r $(CURDIR)/$(LIB_SOURCE_DIR)/hypre-$(hypre_version)/src/hypre/include/ $(CURDIR)/$(LIB_SOURCE_DIR)/hypre

	@# Build PHengLEIv3d0
	@echo "Build PHengLEI CFD executable program..."
	mkdir -p build
	cd build && cmake .. -DCMAKE_INSTALL_PREFIX=$(CURDIR)/build -DCMAKE_CXX_STANDARD_LIBRARIES="-ldl" $(CMAKE_INTEL_FLAG) $(PHengLEI_CMAKE_OPTIONS)
	cd build && make -j$(nProcessor)
endef

# --------------------------- Targets ------------------------------
.PHONY: config distclean

ifdef hpc_mode
config: distclean
	@echo "[HPC Mode] Skipping environment checks"
	$(run-config)
else
config: distclean
	$(check-env)
	$(run-config)
endif

distclean:
	@echo "Cleaning previous build artifacts..."
	rm -rf $(LIB_SOURCE_DIR)/hdf5-$(hdf5_version)
	rm -rf $(LIB_SOURCE_DIR)/CGNS-$(cgns_version)
	rm -rf $(LIB_SOURCE_DIR)/parmetis-$(parmetis_version)
	rm -rf $(LIB_SOURCE_DIR)/metis-$(metis_version)
	rm -rf $(LIB_SOURCE_DIR)/mgrid-$(mgird_version)
	rm -rf $(LIB_SOURCE_DIR)/hypre-$(hypre_version)
	rm -rf build
	@# rm -rf $(LIB_BINARY_DIR)
