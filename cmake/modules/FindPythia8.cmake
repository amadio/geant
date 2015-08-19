# Try to find Pythia8

find_path(PYTHIA8_INCLUDE_DIR Pythia8/Pythia.h PATH_SUFFIXES Pythia8)

find_library(LHAPDFDUMMY_LIBRARY NAMES liblhapdfdummy lhapdfdummy)
find_library(PYTHIA8_LIBRARY NAMES pythia8 libpythia8)

set(PYTHIA8_INCLUDE_DIRS ${PYTHIA8_INCLUDE_DIR})
set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY} ${LHAPDFDUMMY_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia8 DEFAULT_MSG
	LHAPDFDUMMY_LIBRARY PYTHIA8_LIBRARY PYTHIA8_INCLUDE_DIR)

mark_as_advanced(PYTHIA8_INCLUDE_DIR LHAPDFDUMMY_LIBRARY PYTHIA8_LIBRARY)
