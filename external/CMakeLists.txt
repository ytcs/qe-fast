# qe-fast: modified from Quantum ESPRESSO 7.6 (tag qe-7.6), 2026-11-12; see NOTICE.
###########################################################
# FoX
###########################################################
# qe-fast: the FoX submodule declares cmake_minimum_required < 3.5, which CMake 4 refuses
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()
include(fox.cmake)

###########################################################
# WANNIER90
###########################################################
include(wannier90.cmake)

###########################################################
# MBD
###########################################################
include(mbd.cmake)

###########################################################
# DeviceXlib
###########################################################
include(devxlib.cmake)

###########################################################
# Environ
###########################################################
include(environ.cmake)

###########################################################
# D3Q
###########################################################
if("d3q" IN_LIST QE_ENABLE_PLUGINS)
  message(STATUS "Installing D3Q via submodule")
  qe_git_submodule_update(external/d3q)
  add_subdirectory(d3q)
  # qe-fast: externals that d3q calls and QE 7.6 no longer provides (matches, imatches, sgam_ph)
  target_sources(qe_d3q PRIVATE ${PROJECT_SOURCE_DIR}/qefast/compat/d3q_compat.f90)
endif()

###########################################################
# QE-GIPAW
###########################################################
if("gipaw" IN_LIST QE_ENABLE_PLUGINS)
  if(QE_ENABLE_FOX)
    message(STATUS "Installing qe-gipaw via submodule")
    qe_git_submodule_update(external/qe-gipaw)
    add_subdirectory(qe-gipaw)
  else()
    message(FATAL_ERROR "QE-GIPAW requires -DQE_ENABLE_FOX=ON .")
  endif()
endif()
