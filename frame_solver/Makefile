LN_LIB := libfram_solver.so
SHARED_LIBS := libfram_solver.so.0
#SHARED_LIBS := libwrapper_solver.so

TARGETS := $(SHARED_LIBS)
all: $(TARGETS)

#0 hip as backend; 1 cuda as backend;
USE_CUDA := 0
ifeq ($(USE_CUDA), 1)
SRC_BACKEND_DIR := src/nvidia_detail
else
SRC_BACKEND_DIR := src/amd_detail
endif

DEBUG_FLAG = -g
#DEBUG_FLAG = -O3 

SRC_CXX := \
$(SRC_BACKEND_DIR)/hipsolver.cpp \
$(SRC_BACKEND_DIR)/hipsolver_compat.cpp \
$(SRC_BACKEND_DIR)/hipsolver_conversions.cpp \
$(SRC_BACKEND_DIR)/hipsolver_refactor.cpp \
src/common/hipsolver_compat_common.cpp

OBJ_CXX := $(SRC_CXX:.cpp=.o)

FC := gfortran
CXX_FLAGS := $(DEBUG_FLAG) -fPIC -Wall -fopenmp -std=c++17 -D__HIP_PLATFORM_AMD__ 
INC := -I./include -I./include/internal -I./src/include -I/opt/rocm/include
LD_FLAGS := -L/opt/rocm/lib  -lamdhip64 -lrocblas -L../lib -lsolver_syevj
#CXX := hipcc
CXX := g++

INSTALL_PREFIX := /opt/rocm

%.o: %.cpp
	$(CXX)  $(CXX_FLAGS) $(INC) $< -c -o $@

$(SHARED_LIBS): $(OBJ_CXX)
	$(CXX) -g -shared  $^ -o $@ $(LD_FLAGS)
	ln -s $@ $(LN_LIB)


.PHONY: clean
clean:
	rm -rf $(OBJ_CXX) $(TARGETS) $(LN_LIB)


