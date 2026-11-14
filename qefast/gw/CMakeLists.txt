# Copyright (C) 2026 The qe-fast authors
# SPDX-License-Identifier: GPL-2.0-or-later
#
# qe_gw: the GW stages (pair densities, chi0 by Kramers-Kronig, eps^-1, Im W), host and device; chi0.x.
# Application programs outside this tree link qe_gw through -DQEF_EXTRA_DIRS (see the top-level CMakeLists.txt).

set(src_gw ${CMAKE_CURRENT_SOURCE_DIR}/gw_common.f90)
if(QE_ENABLE_CUDA)
    list(APPEND src_gw ${CMAKE_CURRENT_SOURCE_DIR}/chi0_dev.f90)   # cuFFT, cuBLAS
endif()
qe_enable_cuda_fortran("${src_gw}")
qe_add_library(qe_gw ${src_gw})
target_link_libraries(qe_gw
    PRIVATE
        qe_pw
        qe_modules
        qe_upflib
        qe_fftx)
if(QE_ENABLE_CUDA)
    target_link_libraries(qe_gw PRIVATE CUDA::cufft CUDA::cublas)
endif()

# chi0.x: inverse dielectric matrix eps^-1_GG'(q, w) by sum over bands (RPA, local fields)
set(src_chi0_x ${CMAKE_CURRENT_SOURCE_DIR}/chi0.f90)
qe_enable_cuda_fortran("${src_chi0_x}")
qe_add_executable(qe_gw_chi0_exe ${src_chi0_x})
set_target_properties(qe_gw_chi0_exe PROPERTIES OUTPUT_NAME chi0.x)
target_link_libraries(qe_gw_chi0_exe
    PRIVATE
        qe_gw
        qe_pw
        qe_modules
        qe_upflib
        qe_pp
        qe_fftx
        qe_lapack)
if(QE_ENABLE_CUDA)
    target_link_libraries(qe_gw_chi0_exe PRIVATE CUDA::cufft CUDA::cublas)
endif()

qe_install_targets(qe_gw qe_gw_chi0_exe)
