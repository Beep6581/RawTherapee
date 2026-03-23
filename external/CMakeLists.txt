add_library(cjson cJSON/src/cJSON.c)
add_library(cjson::cjson ALIAS cjson)
target_include_directories(cjson
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/cJSON/include
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cJSON/include/cjson
)
