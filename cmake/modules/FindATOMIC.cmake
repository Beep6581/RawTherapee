include(CheckCXXSourceCompiles)

set(ATOMIC_CXX_TEST_SOURCE
"
#include <atomic>
int main() {
	std::atomic<bool> b(false);
	std::atomic<int> i(0);
	std::atomic<unsigned int> ui(0);
	b.exchange(true);
	++i;
	return ++ui;
}
")

set(SAFE_CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")
list(REMOVE_ITEM CMAKE_REQUIRED_LIBRARIES "atomic")

check_cxx_source_compiles("${ATOMIC_CXX_TEST_SOURCE}" HAVE_CXX_ATOMICS_WITHOUT_LIB)

if(NOT HAVE_CXX_ATOMICS_WITHOUT_LIB)
	list(APPEND CMAKE_REQUIRED_LIBRARIES "atomic")
	check_cxx_source_compiles("${ATOMIC_CXX_TEST_SOURCE}" HAVE_CXX_ATOMICS_WITH_LIB)
	if (NOT HAVE_CXX_ATOMICS_WITH_LIB)
		message(FATAL_ERROR "libatomic support needed for std::atomic<> but not 
found.")
	endif()
	set(ATOMIC_LIBRARIES "atomic")
	set(ATOMIC_FOUND ON)
endif()

set(CMAKE_REQUIRED_LIBRARIES "${SAFE_CMAKE_REQUIRED_LIBRARIES}")
unset(SAFE_CMAKE_REQUIRED_LIBRARIES)

mark_as_advanced(ATOMIC_LIBRARIES)
mark_as_advanced(ATOMIC_FOUND)
