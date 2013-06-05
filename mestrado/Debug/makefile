################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf 

# All of the sources participating in the build are defined here
-include sources.mk
-include src/util/subdir.mk
-include src/resolution/grasp/subdir.mk
-include src/resolution/subdir.mk
-include src/problem/subdir.mk
-include src/graph/subdir.mk
-include src/controller/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: MestradoMario

# Tool invocations
MestradoMario: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -pthread -m64 -O2 -fPIC -Wl,-z,noexecstack -Wl,-rpath -Wl,/usr/lib64/mpich2/lib -L/usr/local/lib/boost -L/usr/local/lib -L/usr/lib64/mpich2/lib -o "MestradoMario" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) MestradoMario
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
