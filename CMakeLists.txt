include(${CMAKE_SOURCE_DIR}/cmake/pr_common.cmake)

set(PROJ_NAME util_bsp)
pr_add_library(${PROJ_NAME} STATIC)

pr_add_dependency(${PROJ_NAME} vfilesystem TARGET PUBLIC)
pr_add_dependency(${PROJ_NAME} util_vmf TARGET)
pr_add_dependency(${PROJ_NAME} util_zip TARGET)
pr_add_dependency(${PROJ_NAME} VTFLib13 TARGET)

pr_init_module(${PROJ_NAME})

add_subdirectory(third_party_libs/liblzma)

pr_add_dependency(${PROJ_NAME} lzma TARGET PUBLIC)

pr_finalize(${PROJ_NAME})
