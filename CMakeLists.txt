# cmake version to be used
cmake_minimum_required( VERSION 3.0 )

# project name
project(DifferentiableImage  VERSION 0.0.1 DESCRIPTION "DifferentiableImage library.")

# external libs
find_package(Ceres REQUIRED)
# Or find_package(Ceres REQUIRED PATHS "/path/to/install/dir") if the install directory is not in your PATH

# external libs
# libPeR
find_package(PER REQUIRED per_core per_acquisition_representation per_io per_features per_sensor_pose_estimation per_estimation)
if(PER_FOUND)
	include(${PER_USE_FILE})
endif(PER_FOUND)

# ViSP
find_package(VISP REQUIRED  visp_core visp_robot visp_visual_features visp_io visp_gui visp_vs visp_vision)
if(VISP_FOUND)
	include(${VISP_USE_FILE})
endif(VISP_FOUND)

# Create library 
add_library(${PROJECT_NAME} SHARED src/DifferentiableImage.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC ${CERES_LIBRARIES})

target_include_directories(${PROJECT_NAME} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
    ${CERES_INCLUDE_DIRS}
)


# Create executable
add_executable(DifferentiableImageStandalone src/DifferentiableImageStandalone.cpp)

target_link_libraries(DifferentiableImageStandalone
  PRIVATE
    ${PROJECT_NAME}
)

# Installation 
configure_file("${PROJECT_SOURCE_DIR}/config.cmake.in" "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" @ONLY)
install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake" DESTINATION "share/cmake/${PROJECT_NAME}/")

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

install(DIRECTORY include/ DESTINATION include)

export(TARGETS ${PROJECT_NAME}
    FILE ${PROJECT_NAME}Targets.cmake
)

install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    DESTINATION lib/cmake/${PROJECT_NAME}
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION lib/cmake/${PROJECT_NAME}
)


include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
set(VERSION_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(PROJECT_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake")
set(CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
  ${VERSION_CONFIG}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in" "${PROJECT_CONFIG}"
  INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}"
)

install(FILES "${PROJECT_CONFIG}" "${VERSION_CONFIG}" DESTINATION "${CONFIG_INSTALL_DIR}")

