/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include <casadi/core/casadi_meta.hpp>

namespace casadi {
  const std::string CasadiMeta::version = "3.1.0+";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::git_revision = "0eab93e951b709c091c77772078757dcc0c3d96b";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::git_describe = "3.1.0+139.0eab93e95";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::feature_list = "\n * dynamic-loading , Compile with support for dynamic loading of generated functions (needed for ExternalFunction)\n * hsl-interface , Interface to HSL.\n * lapack-interface , Interface to LAPACK.\n * sundials-interface , Interface to the ODE/DAE integrator suite SUNDIALS.\n * csparse-interface , Interface to the sparse direct linear solver CSparse.\n * blasfeo-interface , Interface to blasfeo.\n * hpmpc-interface , Interface to hpmpc.\n * tinyxml-interface , Interface to the XML parser TinyXML.\n * dsdp-interface , Interface to the interior point SDP solver DSDP.\n * qpoases-interface , Interface to the active-set QP solver qpOASES.\n * blocksqp-interface , Interface to the NLP solver blockSQP.\n * ipopt-interface , Interface to the NLP solver Ipopt.\n";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::build_type = "Release";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler_id = "GNU";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler = "/usr/bin/c++";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::compiler_flags =" -fPIC -fopenmp -DWITH_OPENMP -O3 -DNDEBUG  -DUSE_CXX11 -DHAVE_MKSTEMPS -DCASADI_VERSION=31 -D_USE_MATH_DEFINES -D_SCL_SECURE_NO_WARNINGS -DWITH_DL -DWITH_DEEPBIND";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::modules ="casadi;casadi_linsol_lapacklu;casadi_linsol_lapackqr;casadi_sundials_common;casadi_integrator_cvodes;casadi_integrator_idas;casadi_rootfinder_kinsol;casadi_nlpsol_ipopt;casadi_conic_qpoases;casadi_linsol_csparse;casadi_linsol_csparsecholesky;casadi_linsol_ma27;casadi_xmlfile_tinyxml;casadi_nlpsol_blocksqp;casadi_conic_hpmpc;casadi_conic_nlpsol;casadi_importer_shell;casadi_integrator_rk;casadi_integrator_collocation;casadi_interpolant_linear;casadi_linsol_symbolicqr;casadi_nlpsol_sqpmethod;casadi_nlpsol_scpgen;casadi_rootfinder_newton;casadi_rootfinder_nlpsol";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::plugins ="Linsol::lapacklu;Linsol::lapackqr;Integrator::cvodes;Integrator::idas;Rootfinder::kinsol;Nlpsol::ipopt;Conic::qpoases;Linsol::csparse;Linsol::csparsecholesky;Linsol::ma27;XmlFile::tinyxml;Nlpsol::blocksqp;Conic::hpmpc;Conic::nlpsol;Importer::shell;Integrator::rk;Integrator::collocation;Interpolant::linear;Linsol::symbolicqr;Nlpsol::sqpmethod;Nlpsol::scpgen;Rootfinder::newton;Rootfinder::nlpsol";  // NOLINT(whitespace/line_length)
  const std::string CasadiMeta::install_prefix ="/home/greghorn/casadi_install";  // NOLINT(whitespace/line_length)
} // namespace casadi
