#  James Bigler, NVIDIA Corp (nvidia.com - jbigler)
#
#  Copyright (c) 2008 - 2009 NVIDIA Corporation.  All rights reserved.
#
#  This code is licensed under the MIT License.  See the FindCUDA.cmake script
#  for the text of the license.

# The MIT License
#
# License for the specific language governing rights and limitations under
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


##########################################################################
# This file runs the nvcc commands to produce the desired output file along with
# the dependency file needed by CMake to compute dependencies.  In addition the
# file checks the output of each command and if the command fails it deletes the
# output files.

# Input variables
#
# verbose:BOOL=<>          OFF: Be as quiet as possible (default)
#                          ON : Describe each step
#
# build_configuration:STRING=<> Typically one of Debug, MinSizeRel, Release, or
#                               RelWithDebInfo, but it should match one of the
#                               entries in CUDA_HOST_FLAGS. This is the build
#                               configuration used when compiling the code.  If
#                               blank or unspecified Debug is assumed as this is
#                               what CMake does.
#
# generated_file:STRING=<> File to generate.  This argument must be passed in.
#
# generated_cubin_file:STRING=<> File to generate.  This argument must be passed
#                                                   in if build_cubin is true.

if(NOT generated_host_obj_file)
  message(FATAL_ERROR "You must specify generated_host_obj_file on the command line")

endif()

# Set these up as variables to make reading the generated file easier
set(CMAKE_COMMAND "@CMAKE_COMMAND@") # path
set(object_source_files "@object_source_files@") # path
set(CUDA_make2cmake "@CUDA_make2cmake@") # path
set(CUDA_parse_cubin "@CUDA_parse_cubin@") # path

# We won't actually use these variables for now, but we need to set this, in
# order to force this file to be run again if it changes.
set(generated_file_path "@generated_file_path@") # path
set(generated_file_internal "@generated_host_obj_file@") # path

set(CUDA_NVCC_EXECUTABLE "@CUDA_NVCC_EXECUTABLE@") # path

#message("*************************************")
#message("nvlink")
#message(${object_source_files})
#message("*************************************")


# cuda_execute_process - Executes a command with optional command echo and status message.
#
#   status  - Status message to print if verbose is true
#   command - COMMAND argument from the usual execute_process argument structure
#   ARGN    - Remaining arguments are the command with arguments
#
#   CUDA_result - return value from running the command
#
# Make this a macro instead of a function, so that things like RESULT_VARIABLE
# and other return variables are present after executing the process.
macro(cuda_execute_process status command)
  set(_command ${command})
  if(NOT _command STREQUAL "COMMAND")
    message(FATAL_ERROR "Malformed call to cuda_execute_process.  Missing COMMAND as second argument. (command = ${command})")
  endif()
  if(verbose)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E echo -- ${status})
    # Now we need to build up our command string.  We are accounting for quotes
    # and spaces, anything else is left up to the user to fix if they want to
    # copy and paste a runnable command line.
    set(cuda_execute_process_string)
    foreach(arg ${ARGN})
      # If there are quotes, excape them, so they come through.
      string(REPLACE "\"" "\\\"" arg ${arg})
      # Args with spaces need quotes around them to get them to be parsed as a single argument.
      if(arg MATCHES " ")
        list(APPEND cuda_execute_process_string "\"${arg}\"")
      else()
        list(APPEND cuda_execute_process_string ${arg})
      endif()
    endforeach()
    # Echo the command
    execute_process(COMMAND ${CMAKE_COMMAND} -E echo ${cuda_execute_process_string})
  endif(verbose)
  # Run the command
  execute_process(COMMAND ${ARGN} RESULT_VARIABLE CUDA_result )
endmacro()

# Delete the target file
cuda_execute_process(
  "Removing ${generated_host_obj_file}"
  COMMAND "${CMAKE_COMMAND}" -E remove "${generated_host_obj_file}"
  )

cuda_execute_process(
  "Generating ${generated_host_obj_file}"
  COMMAND "${CUDA_NVCC_EXECUTABLE}"
  -arch=sm_20
  -dlink
  "${object_source_files}"
  -o "${generated_host_obj_file}"
  )


if(CUDA_result)
  # Since nvcc can sometimes leave half done files make sure that we delete the output file.
  cuda_execute_process(
    "Removing ${generated_host_obj_file}"
    COMMAND "${CMAKE_COMMAND}" -E remove "${generated_host_obj_file}"
    )
  message(FATAL_ERROR "Error generating file ${generated_host_obj_file}")
else()
  if(verbose)
    message("Generated ${generated_host_obj_file} successfully.")
  endif()
endif()






