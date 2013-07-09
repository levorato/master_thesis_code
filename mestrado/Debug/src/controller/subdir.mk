################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/controller/CommandLineInterfaceController.cpp \
../src/controller/SimpleTextGraphFileReader.cpp 

OBJS += \
./src/controller/CommandLineInterfaceController.o \
./src/controller/SimpleTextGraphFileReader.o 

CPP_DEPS += \
./src/controller/CommandLineInterfaceController.d \
./src/controller/SimpleTextGraphFileReader.d 


# Each subdirectory must supply rules for building sources it contributes
src/controller/%.o: ../src/controller/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/mpich2-x86_64 -I/usr/lib64/mpich2/lib -I/usr/local/include -I/usr/include/mpi -I/usr/include -I$HOME/include -O2 -g3 -Wall -c -fmessage-length=0 -pthread -m64 -Wl,-z,noexecstack -Wl,-rpath -Wl,/usr/lib64/mpich2/lib -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


