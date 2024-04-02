EXE_SYEVJ_DIR 	:= examples/
LIB_SYEVJ_DIR 	:= solver_syevj/
LIB_FRAME_DIR 	:= frame_solver/

LIB_DIR		:= ./lib

EXE_SYEVJ 	:= app_solver_Dsyevj
LIB_SYEVJ	:= libsolver_syevj.so
LIB_FRAME	:= libfram_solver.so

all: $(EXE_SYEVJ)


$(EXE_SYEVJ):  $(LIB_SYEVJ) $(LIB_FRAME)
	$(MAKE) -j -C $(EXE_SYEVJ_DIR)
	cp $(EXE_SYEVJ_DIR)/$(EXE_SYEVJ) $(PWD)/

$(LIB_FRAME): $(LIB_SYEVJ)
	$(MAKE) -j -C $(LIB_FRAME_DIR)
	mkdir -p $(LIB_DIR)
	cp $(LIB_FRAME_DIR)/$(LIB_FRAME) $(LIB_DIR)/

$(LIB_SYEVJ):
	$(MAKE) -j -C $(LIB_SYEVJ_DIR)
	mkdir -p $(LIB_DIR)
	cp $(LIB_SYEVJ_DIR)/$(LIB_SYEVJ) $(LIB_DIR)/



.PHONY: clean
clean:
	rm -rf $(EXE_SYEVJ)
	rm -rf $(LIB_DIR)
	$(MAKE) clean -C $(EXE_SYEVJ_DIR)
	$(MAKE) clean -C $(LIB_FRAME_DIR)
	$(MAKE) clean -C $(LIB_SYEVJ_DIR)
	
.PHONY: echo
	@echo "Hello"