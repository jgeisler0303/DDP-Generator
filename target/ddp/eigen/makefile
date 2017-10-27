OPTDEF_NAME:= $(PROBLEM_FILE:%.mac=%)
PROBLEM_NAME:= $(OPTDEF_NAME:optDef%=%)
PROGRAM:= $(PROBLEM_NAME:%=ddp%_eigen_cli)

PROBLEM_FILE:=  optDef$(PROBLEM_NAME).mac

GEN_DIR:= gendir_$(PROBLEM_NAME)_ddp_eigen
TEST_DIR:= testdir_$(PROBLEM_NAME)_ddp_eigen
WORKING_DIR := ${CURDIR}

RESULT_FILE:=  $(TEST_DIR)/result_$(PROBLEM_NAME).dat
RESULT_FILE_:=  $(TEST_DIR)/result_%.dat
RESULT_EPS:=  $(TEST_DIR)/result_$(PROBLEM_NAME).eps
ITERATION_FILE:=  $(TEST_DIR)/iterations_$(PROBLEM_NAME).dat
ITERATION_FILE_:=  $(TEST_DIR)/iterations_%.dat
ITERATION_EPS:=  $(TEST_DIR)/iterations_$(PROBLEM_NAME).eps

USER_FLAGS_MK:= $(OPTDEF_NAME)_eigen.mk

ifeq ($(wildcard $(PROBLEM_FILE)),)
$(error Problem file $(PROBLEM_FILE) not found)
endif

-include $(USER_FLAGS_MK)

OBJ_DIR:= $(GEN_DIR)/obj

MKDIR_P:= mkdir -p
CPP:= g++
INCLUDES:= -I$(DDP_BASE_DIR) -I$(GEN_DIR) -I$(DDP_BASE_DIR)/../../../util/inih -I$(DDP_BASE_DIR)/../../../util
CPPFLAGS:= $(USER_FLAGS) -Wall -std=c++0x 
LIBS= -lm

TARGET_OBJS:= back_pass.o boxQP.o iterate.o line_search.o printMat.o ddp_cli.o
GEN_OBJS:= iLQG_func.o
CLI_OBJS:= strtonum.o ini.o
OBJS_:= $(TARGET_OBJS) $(GEN_OBJS) $(CLI_OBJS) 
OBJS:= $(addprefix $(GEN_DIR)/obj/, $(OBJS_))

all : cli

cli : $(PROGRAM)

test : $(RESULT_EPS) $(ITERATION_EPS)

# suppress makes desire to remake these includes
$(OBJS:.o=.d) : ;
# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

define COMPILE_CMD=
$(CPP) -MMD -MT $@ -MP -c $(CPPFLAGS) $(INCLUDES) $< -o $@
endef

$(PROGRAM) : $(OBJS)
	$(CPP) $(OBJS) -o $(PROGRAM) $(LIBS)

$(GEN_DIR)/obj/%.o : $(GEN_DIR)/%.cpp $(GEN_DIR)/iLQG_problem.h
	$(COMPILE_CMD)

$(GEN_DIR)/obj/%.o : $(DDP_BASE_DIR)/%.cpp $(GEN_DIR)/iLQG_problem.h
	$(COMPILE_CMD)

$(GEN_DIR)/obj/%.o : $(DDP_BASE_DIR)/../../../util/%.c
	$(COMPILE_CMD)

$(GEN_DIR)/obj/%.o : $(DDP_BASE_DIR)/../../../util/inih/%.c
	$(COMPILE_CMD)


$(GEN_DIR)/%_func.cpp $(GEN_DIR)/%_problem.h : $(PROBLEM_FILE) | $(GEN_DIR)
	@$(RM) $(GEN_DIR)/maxima_complete; \
	cd $(DDP_BASE_DIR); \
	maxima -q --batch-string="load(\"make_iLQG.mac\")\$$ gen_ddp(\"$(WORKING_DIR)/$(PROBLEM_FILE)\", \"$(WORKING_DIR)/$(GEN_DIR)\")\$$" 2>$(WORKING_DIR)/$(GEN_DIR)/maxima_errors.log && \
	if [ -s  "$(WORKING_DIR)/$(GEN_DIR)/maxima_errors.log" ] || ! [ -e  "$(WORKING_DIR)/$(GEN_DIR)/maxima_complete" ]; \
	then mv -f $(WORKING_DIR)/$(GEN_DIR)/iLQG_func.cpp $(WORKING_DIR)/$(GEN_DIR)/faulty_iLQG_func.cpp; \
	mv -f $(WORKING_DIR)/$(GEN_DIR)/iLQG_problem.h $(WORKING_DIR)/$(GEN_DIR)/faulty_iLQG_problem.h; \
	$(RM) $(GEN_DIR)/maxima_complete; \
	exit 1; fi;

clean : clean_cli

clean_cli :
	$(RM) $(PROGRAM)
	$(RM) -r $(GEN_DIR)/obj

$(OBJS) : $(USER_FLAGS_MK) | $(OBJ_DIR)

$(OBJ_DIR) : | $(GEN_DIR)
	$(MKDIR_P) $(OBJ_DIR)

$(GEN_DIR) :
	$(MKDIR_P) $(GEN_DIR)

$(TEST_DIR) :
	$(MKDIR_P) $(TEST_DIR)

$(RESULT_FILE_) $(ITERATION_FILE_) : $(PROGRAM) | $(TEST_DIR)
	./$(PROGRAM) > /dev/null
	mv ddp_result.dat $(RESULT_FILE)
	mv ddp_iter.dat $(ITERATION_FILE)

$(RESULT_EPS) : $(RESULT_FILE)
	gnuplot -e "outfile='$(RESULT_EPS)'" -e "datafile='$(RESULT_FILE)'" ddp$(PROBLEM_NAME).plt

$(ITERATION_EPS) : $(ITERATION_FILE)
	gnuplot -e "outfile='$(ITERATION_EPS)'" -e "datafile='$(ITERATION_FILE)'" $(DDP_BASE_DIR)/../common/ddpConvergence.plt

$(PROBLEM_FILE) : ;

$(USER_FLAGS_MK) : 
	@echo \# Here you can specify additional compiler flags > $(USER_FLAGS_MK)
	@echo \#USER_FLAGS:= -O2 -g -march=native -mfpu=neon -mfloat-abi=hard -DEIGEN_NO_MALLOC -DEIGEN_NO_DEBUG -DEIGEN_DONT_PARALLELIZE >> $(USER_FLAGS_MK)

.PHONY : clean all cli clean_cli 

.SECONDARY : $(GEN_DIR)/iLQG_func.c $(GEN_DIR)/iLQG_problem.h