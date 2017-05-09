# Module for locating libhwloc
#
# Read-only variables:
#   hwloc_FOUND
#     Indicates that the library has been found.
#
#   HWLOC_LIBRARIES
#     Points to the d
#     The content of this variable can be passed to target_link_libraries.
#

# Use pkg-config to fetch the contents of the .pc file
if(NOT PKGCONFIG_FOUND)
  find_package(PkgConfig QUIET)
endif()

if(HWLOC_ROOT)
  set(ENV{PKG_CONFIG_PATH} "${HWLOC_ROOT}/lib/pkgconfig")
else()
  foreach(PREFIX ${CMAKE_PREFIX_PATH})
    set(PKG_CONFIG_PATH "${PKG_CONFIG_PATH}:${PREFIX}/lib/pkgconfig")
  endforeach()
  set(ENV{PKG_CONFIG_PATH} "${PKG_CONFIG_PATH}:$ENV{PKG_CONFIG_PATH}")
endif()

if(hwloc_FIND_REQUIRED)
  set(_hwloc_OPTS "REQUIRED")
elseif(hwloc_FIND_QUIETLY)
  set(_hwloc_OPTS "QUIET")
else()
  set(_hwloc_output 1)
endif()

if(hwloc_FIND_VERSION)
  if(hwloc_FIND_VERSION_EXACT)
    pkg_check_modules(HWLOC ${_hwloc_OPTS} hwloc=${hwloc_FIND_VERSION})
  else()
    pkg_check_modules(HWLOC ${_hwloc_OPTS} hwloc>=${hwloc_FIND_VERSION})
  endif()
else()
  pkg_check_modules(HWLOC ${_hwloc_OPTS} hwloc)
endif()

if(HWLOC_FOUND)
#  include(FindPackageHandleStandardArgs)
#  find_package_handle_standard_args(hwloc DEFAULT_MSG hwloc_LIBRARIES hwloc_INCLUDE_DIR)
  set(HWLOC_INCLUDE_DIRS ${HWLOC_INCLUDEDIR})
#  set(HWLOC_LIBRARIES ${HWLOC_LIBDIR}/libhwloc)
endif()
