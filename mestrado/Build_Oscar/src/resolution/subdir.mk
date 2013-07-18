################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/resolution/ProblemSolverFacade.cpp \
../src/resolution/ResolutionStrategy.cpp 

OBJS += \
./src/resolution/ProblemSolverFacade.o \
./src/resolution/ResolutionStrategy.o 

CPP_DEPS += \
./src/resolution/ProblemSolverFacade.d \
./src/resolution/ResolutionStrategy.d 


# Each subdirectory must supply rules for building sources it contributes
src/resolution/%.o: ../src/resolution/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	icpc -DBOOST_ALL_DYN_LINK -I/usr/include/mpi -DMPIBULL_IGNORE_CXX_SEEK -DHAVE_MPI_CXX -I/opt/mpi/mpibull2-1.3.9-18.s/include -L/opt/mpi/mpibull2-1.3.9-18.s/lib/drivers/osock -L/opt/mpi/mpibull2-1.3.9-18.s/lib -lmpibinding_cxx -lmpi -lrt -ldl -L/opt/mpi/mpibull2-1.3.9-18.s/lib/pmi -lpmi -lpthread -luuid -I/usr/local/include -I$$HOME/include -O2 -g3 -Wall -c -fmessage-length=0 -pthread -m64 -Wl,-z,noexecstack -Wl,-rpath -Wl,/usr/lib64/mpich2/lib -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


