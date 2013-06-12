################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/problem/CCProblem.cpp \
../src/problem/ClusteringProblem.cpp \
../src/problem/ClusteringProblemFactory.cpp \
../src/problem/RCCProblem.cpp 

OBJS += \
./src/problem/CCProblem.o \
./src/problem/ClusteringProblem.o \
./src/problem/ClusteringProblemFactory.o \
./src/problem/RCCProblem.o 

CPP_DEPS += \
./src/problem/CCProblem.d \
./src/problem/ClusteringProblem.d \
./src/problem/ClusteringProblemFactory.d \
./src/problem/RCCProblem.d 


# Each subdirectory must supply rules for building sources it contributes
src/problem/%.o: ../src/problem/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/mpi -I/usr/lib64/mpich2/lib -I/usr/lib64/mpich2/lib -I/usr/local/include -I$$HOME/include -O2 -g3 -Wall -c -fmessage-length=0 -pthread -m64 -Wl,-z,noexecstack -Wl,-rpath -Wl,/usr/lib64/mpich2/lib -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


