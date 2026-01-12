cmake_minimum_required (VERSION 3.18)

set(CMAKE_CUDA_ARCHITECTURES 89)

if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("MolecularDynamicsEngine" LANGUAGES CXX CUDA)

find_package(CUDAToolkit REQUIRED)

link_directories("C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.1/lib/x64")

add_executable (MolecularDynamicsEngine
	"main.cpp"
	"Box.h"
	"Box.cpp"
	"Box.cu"
	"Matrix.h"
)

target_link_libraries(MolecularDynamicsEngine CUDA::cudart)

set_target_properties(MolecularDynamicsEngine PROPERTIES
	CXX_STANDARD 20
	CUDA_STANDARD 17
)