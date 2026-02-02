cmake_minimum_required (VERSION 3.18)

set(CMAKE_CUDA_ARCHITECTURES 89)
set(CMAKE_CUDA_COMPILER "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.6/bin/nvcc.exe")

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler")

project ("MolecularDynamicsEngine" LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)

add_executable (MolecularDynamicsEngine
	"main.cpp"
	"Box.h"
	"Box.cpp"
	"Box.cu"
)

target_link_libraries(MolecularDynamicsEngine CUDA::cudart)

set_target_properties(MolecularDynamicsEngine PROPERTIES
	CXX_STANDARD 20
	CUDA_STANDARD 17
)
