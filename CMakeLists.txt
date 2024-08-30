cmake_minimum_required(VERSION 3.5)

project(PassivityTorqueFeedback LANGUAGES CXX VERSION 1.0.0)

include(CTest)
enable_testing()

# Find mc_rtc if we are not building inside mc_rtc itself
if(NOT TARGET mc_rtc::mc_control)
  find_package(mc_rtc REQUIRED)
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

# Configure and install the plugin configuration
configure_file(etc/PassivityTorqueFeedback.in.yaml "${CMAKE_CURRENT_BINARY_DIR}/etc/PassivityTorqueFeedback.yaml")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/etc/PassivityTorqueFeedback.yaml" DESTINATION "${MC_PLUGINS_RUNTIME_INSTALL_PREFIX}/etc")

set(plugin_SRC
  src/PassivityTorqueFeedback.cpp
  )

set(plugin_HDR
  src/PassivityTorqueFeedback.h
  include/utils/ROSSubscriber.h

)


# The add_plugin macro is provided by mc_rtc
add_plugin(PassivityTorqueFeedback "${plugin_SRC}" "${plugin_HDR}")

find_package(jrl-qp REQUIRED)
target_link_libraries(PassivityTorqueFeedback PUBLIC jrl-qp::jrl-qp)
target_link_libraries(PassivityTorqueFeedback PUBLIC mc_rtc::mc_rtc_ros)


if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
