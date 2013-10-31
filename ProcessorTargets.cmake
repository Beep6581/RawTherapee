# We have to create a label variable if we want to display it in AboutThisBuild.txt...

# This first choice should be used for official releases
set(PROC_TARGET_1_LABEL generic x86 CACHE STRING "Processor-1 label - should be used for official Windows release")
set(PROC_TARGET_1_FLAGS "-mtune=generic" CACHE STRING "Processor-1 flags")

# This second choice should be used for your own build only
set(PROC_TARGET_2_LABEL native CACHE STRING "Processor-2 label - use it for your own build")
set(PROC_TARGET_2_FLAGS "-march=native" CACHE STRING "Processor-2 flags")

# The later choices is intended to be used if you want to provide specific builds, but it should match your own processor
# You can cross compile but you have to know what you're doing, this mechanism has not been designed for that

set(PROC_TARGET_3_LABEL pentium CACHE STRING "Processor-3 label - use it to provide a pentium optimized build, if you have this processor")
set(PROC_TARGET_3_FLAGS "-march=pentium" CACHE STRING "Processor-3 flags")

set(PROC_TARGET_4_LABEL pentium4 CACHE STRING "Processor-4 label - use it to provide a pentium4 optimized build, if you have this processor")
set(PROC_TARGET_4_FLAGS "-march=pentium4" CACHE STRING "Processor-4 flags")

set(PROC_TARGET_5_LABEL core2 CACHE STRING "Processor-5 label - use it to provide a core2 optimized build, if you have this processor")
set(PROC_TARGET_5_FLAGS "-march=core2" CACHE STRING "Processor-5 flags")

set(PROC_TARGET_6_LABEL corei7 CACHE STRING "Processor-6 label - use it to provide a corei7 optimized build, if you have this processor")
set(PROC_TARGET_6_FLAGS "-march=corei7" CACHE STRING "Processor-6 flags")

set(PROC_TARGET_7_LABEL athlon-4 CACHE STRING "Processor-7 label - use it to provide a athlon-4 optimized build, if you have this processor")
set(PROC_TARGET_7_FLAGS "-march=athlon-4" CACHE STRING "Processor-7 flags")

set(PROC_TARGET_8_LABEL athlon64 CACHE STRING "Processor-8 label - use it to provide a athlon64 optimized build, if you have this processor")
set(PROC_TARGET_8_FLAGS "-march=athlon64" CACHE STRING "Processor-8 flags")

set(PROC_TARGET_9_LABEL phenomX4 CACHE STRING "Processor-9 label - use it to provide a phenomX4 optimized build, if you have this processor")
set(PROC_TARGET_9_FLAGS "-march=amdfam10" CACHE STRING "Processor-9 flags")

#set(PROC_TARGET__LABEL procLabel CACHE STRING "Processor- label")
#set(PROC_TARGET__FLAGS "procFlags" CACHE STRING "Processor- flags")
