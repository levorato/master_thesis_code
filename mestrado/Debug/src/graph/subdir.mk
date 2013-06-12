################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/graph/Clustering.cpp \
../src/graph/Graph.cpp \
../src/graph/SequentialNeighborhoodGen.cpp 

OBJS += \
./src/graph/Clustering.o \
./src/graph/Graph.o \
./src/graph/SequentialNeighborhoodGen.o 

CPP_DEPS += \
./src/graph/Clustering.d \
./src/graph/Graph.d \
./src/graph/SequentialNeighborhoodGen.d 


# Each subdirectory must supply rules for building sources it contributes
src/graph/%.o: ../src/graph/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/mpi -I/usr/lib64/mpich2/lib -I/usr/lib64/mpich2/lib -I/usr/local/include -I$$HOME/include -O2 -g3 -Wall -c -fmessage-length=0 -pthread -m64 -Wl,-z,noexecstack -Wl,-rpath -Wl,/usr/lib64/mpich2/lib -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


