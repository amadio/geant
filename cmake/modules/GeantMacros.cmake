MACRO(add_headers headers)
  set(allh)
  foreach(_header ${headers})		     
    GET_FILENAME_COMPONENT(_base_header ${_header} NAME)
    add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/include/${_base_header}
    COMMAND cp ${_header} ${PROJECT_SOURCE_DIR}/include/${_base_header}
    DEPENDS ${_header})
    set(allh ${allh} ${PROJECT_SOURCE_DIR}/include/${_base_header})
  endforeach(_header)
  add_custom_command(OUTPUT ${PROJECT_SOURCE_DIR}/include
  		   COMMAND mkdir ${PROJECT_SOURCE_DIR}/include)
  string(REPLACE ${PROJECT_SOURCE_DIR}/ "" _leaf ${CMAKE_CURRENT_SOURCE_DIR})
  add_custom_target(${_leaf}_headers DEPENDS ${PROJECT_SOURCE_DIR}/include ${allh})
  add_dependencies(_headers ${_leaf}_headers)
endmacro(add_headers)

