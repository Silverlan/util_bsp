include(${CMAKE_SOURCE_DIR}/cmake/pr_common.cmake)

set(PROJ_NAME util_bsp)
pr_add_library(${PROJ_NAME} STATIC)

pr_add_dependency(${PROJ_NAME} vfilesystem TARGET PUBLIC)
pr_add_dependency(${PROJ_NAME} util_vmf TARGET)
pr_add_dependency(${PROJ_NAME} util_zip TARGET)
pr_add_external_dependency(${PROJ_NAME} vtflib_build HEADER_ONLY)
pr_add_external_dependency(${PROJ_NAME} VTFLib LIBRARY)

pr_add_headers(${PROJ_NAME} "include/")
pr_add_sources(${PROJ_NAME} "src/")

add_subdirectory(third_party_libs/liblzma)

pr_add_dependency(${PROJ_NAME} lzma TARGET PUBLIC)

pr_finalize(${PROJ_NAME})
