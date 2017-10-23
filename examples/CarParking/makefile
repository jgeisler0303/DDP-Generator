DDP_BASE_DIR:= /media/fritzbox/home/jgeisler/git/DDP-Generator
TARGETDIR:= $(DDP_BASE_DIR)/target/ddp/eigen
GEN_DIR:= gendir_Car_ddp_eigen
PROGRAM:= ddpCar_cli


CPP:= g++
INCLUDES:= -I$(TARGETDIR) -I$(GEN_DIR) -I$(DDP_BASE_DIR)/inih -I$(DDP_BASE_DIR)/util
CPPFLAGS:= -O2 -std=c++0x -ggdb

TARGET_OBJS:= back_pass.o boxQP.o iterate.o line_search.o printMat.o ddp_cli.o
GEN_OBJS:= iLQG_func.o
CLI_OBJS:= strtonum.o ini.o
OBJS_:= $(TARGET_OBJS) $(GEN_OBJS) $(CLI_OBJS) 
OBJS:= $(addprefix $(GEN_DIR)/obj/, $(OBJS_))


all : cli

cli : $(OBJS)
	$(CPP) $(OBJS) -o $(PROGRAM)

$(GEN_DIR)/obj/%.o : $(GEN_DIR)/%.cpp
	$(CPP) -c $(CPPFLAGS) $(INCLUDES) $< -o $@


$(GEN_DIR)/obj/%.o : $(TARGETDIR)/%.cpp
	$(CPP) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

$(GEN_DIR)/obj/%.o : $(DDP_BASE_DIR)/util/%.c
	$(CPP) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

$(GEN_DIR)/obj/%.o : $(DDP_BASE_DIR)/inih/%.c
	$(CPP) -c $(CPPFLAGS) $(INCLUDES) $< -o $@

.PHONY: clean all
