add_library(cjson cJSON/src/cJSON.c)
add_library(cjson::cjson ALIAS cjson)
target_include_directories(cjson
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/cJSON/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cJSON/include/cjson
)
rt_setup_target(cjson)

if(NOT WITH_SYSTEM_KLT)
    add_library(klt
        klt/src/convolve.cc
        klt/src/error.cc
        klt/src/klt_util.cc
        klt/src/klt.cc
        klt/src/pnmio.cc
        klt/src/pyramid.cc
        klt/src/selectGoodFeatures.cc
        klt/src/storeFeatures.cc
        klt/src/trackFeatures.cc
        klt/src/writeFeatures.cc
    )
    add_library(klt::klt ALIAS klt)
    target_include_directories(klt
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/klt/include
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/klt/include/klt
    )
    rt_setup_target(klt)
else()
    add_library(klt INTERFACE)
    add_library(klt::klt ALIAS klt)
    target_link_libraries(klt INTERFACE PkgConfig::KLT)
endif()
