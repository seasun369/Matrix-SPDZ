cmake_minimum_required(VERSION 3.16)
project(VOLE)

add_executable(main main.cpp)
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_CXX_FLAGS "-O3")



find_package(libOTe REQUIRED HINTS /user/local
	COMPONENTS
		sodium
		boost
		openssl
		#circuits
		#sse
		#avx
		#no_sse
		#no_avx
		bitpolymul
		simplestot
		#mrr
		#mrr_twist
		#mr
		mr_kyber
		kos
		iknp
		silentot
		softspoken_ot
		delta_kos
		silent_vole
		oos
		kkrt
)
target_compile_options(main PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-std=c++20>)
target_link_libraries(main oc::libOTe oc::libOTe_Tests oc::cryptoTools)

