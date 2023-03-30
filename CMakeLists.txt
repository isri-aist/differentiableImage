# cmake version to be used
cmake_minimum_required( VERSION 3.0 )

# project name
project(differentiableImage)

# external libs
find_package(Ceres REQUIRED)
# Or find_package(Ceres REQUIRED PATHS "/path/to/install/dir") if the install directory is not in your PATH

# external libs
find_package(OpenCV REQUIRED)

# target multiSalign_sphericalHomog_hyperspectral2gray
add_executable(differentiableImage differentiableImage.cpp )

SET(OpenCV_INCLUDE_SOURCE_DIRS "${OpenCV_INCLUDE_DIRS}/opencv2/"
																"${OpenCV_SRC_DIR}/modules/calib3d/src/"
																"${OpenCV_SRC_DIR}/modules/core/include/")

target_include_directories(differentiableImage
  PRIVATE
    ${CERES_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_SOURCE_DIRS}
    ${libnpy_INCLUDE_DIRS}
)

target_link_libraries(differentiableImage
  PRIVATE
    ${CERES_LIBRARIES}
    ${OpenCV_LIBS}
)

