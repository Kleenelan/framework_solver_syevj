SRC := app_solver_Dsyevj.cpp
EXE := app_solver_Dsyevj

all: $(EXE)

#LD_SOLVER := -L../lib -lfram_solver -lsolver_syevj
#LD_ROCBLAS := -L/home/hipper/ex_hipsolver_rocsolver_suitesparse/local_d/rocblas/lib -lrocblas
#INC      := -I/opt/rocm/include
#LD_FLAGS := -L/opt/rocm/lib -lamdhip64 $(LD_ROCBLAS) $(LD_SOLVER) -D__HIP_PLATFORM_AMD__  #-lrocsolver

INC_CUDA := -I /usr/local/cuda/include
LD_CUDA  := -L /usr/local/cuda/lib64 -lcudart -lcublas -lcusolver

%: %.cpp
	g++ -g $< $(INC_CUDA) -o $@ $(LD_CUDA)

#	g++ -g $< -o $@ $(INC) $(LD_FLAGS)


.PHONY: clean
clean:
	rm -rf $(EXE)



