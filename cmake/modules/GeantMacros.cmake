 # define a function downloading from a URL into a local file LOCALFILE
 function(FILE_DOWNLOAD FILE_URL LOCALFILE )
  if(APPLE)
      #execute_process(COMMAND curl ${FILE_URL} -o  ${LOCALFILE})
      file(DOWNLOAD ${FILE_URL} ${LOCALFILE})
  else()
     execute_process(COMMAND wget -q ${FILE_URL} -O ${LOCALFILE})
     #file(DOWNLOAD ${FILE_URL} ${FILE_URL})
  endif()
 endfunction(FILE_DOWNLOAD)
 # end of function FILE DOWNLOAD

 # define a function checking md5 hashes
 # result is stored in MD5MATCHES ( 1 == true, 0 == false )
 function(CHECKMD5 FILETOCHECK EXPECTEDMD5HASH MD5MATCHES)
     if(APPLE)
         execute_process(COMMAND md5 ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
         string(LENGTH ${MD5SUM} MD5LENGTH)
         MATH(EXPR START "${MD5LENGTH} - 33")
         string(SUBSTRING ${MD5SUM} ${START} 32 MD5SUM)
     else()
         execute_process(COMMAND md5sum ${FILETOCHECK} OUTPUT_VARIABLE MD5SUM)
         string(SUBSTRING ${MD5SUM} 0 32 MD5SUM)
     endif()
     if(MD5SUM STREQUAL EXPECTEDMD5HASH)
       set(${MD5MATCHES} 1 PARENT_SCOPE)
     else()
       set(${MD5MATCHES} 0 PARENT_SCOPE)
     endif()
 endfunction(CHECKMD5)

 # actual function for managing the download
 function(DOWNLOAD_IF_NOT_INSTALLED FILE_URL LOCALFILE TARGETPATH MD5HASH )
   find_file(FOUNDFILE NAMES ${LOCALFILE} PATHS ${TARGETPATH})
   if(FOUNDFILE STREQUAL "FOUNDFILE-NOTFOUND")
       # set need download
       message(STATUS "need download of ${LOCALFILE} since not found")
       set( NEEDTODOWNLOAD 1 )
   else()
       # check md5
       message(STATUS "found existing file ${LOCALFILE}")
       CHECKMD5( ${FOUNDFILE} ${MD5HASH} MD5CORRECT )
       if( ${MD5CORRECT} STREQUAL "1" )
           # do not set download flag
           set( NEEDTODOWNLOAD 0 )
       else( )
           # set need download
           message(STATUS "hash ${MD5HASH} not correct for file ${FOUNDFILE} ${MD5CORRECT}" )
           set( NEEDTODOWNLOAD 1 )
       endif( )
   endif()

   if( ${NEEDTODOWNLOAD} STREQUAL 1 )
       message(STATUS " downloading ... ")
       set(DOWNLOADLOCATION "${TARGETPATH}/${LOCALFILE}")
       FILE_DOWNLOAD( ${FILE_URL} ${DOWNLOADLOCATION} )
   else()
       message(STATUS " doing nothing ... ")
   endif()
   # in principle have to check now if download succeeded and has right MD5
   # TOBEDONE

   # this is annoying but we have to clear FOUNDFILE SINCE THIS IS TREATED LIKE A STATIC VARIABLE
   unset(FOUNDFILE CACHE)
 endfunction(DOWNLOAD_IF_NOT_INSTALLED)



