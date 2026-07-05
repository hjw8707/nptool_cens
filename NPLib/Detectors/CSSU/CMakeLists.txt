add_custom_command(OUTPUT TCSSUDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCSSUData.h TCSSUDataDict.cxx TCSSUData.rootmap libNPCSSU.dylib DEPENDS TCSSUData.h)
add_library(NPCSSU SHARED TCSSUData.cxx TCSSUDataDict.cxx)
target_link_libraries(NPCSSU ${ROOT_LIBRARIES} NPCore)
install(FILES TCSSUData.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
