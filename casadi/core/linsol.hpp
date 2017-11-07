/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_LINSOL_HPP
#define CASADI_LINSOL_HPP

#include "function.hpp"
#include "printable.hpp"

namespace casadi {

  // Forward declaration of internal class
  class LinsolInternal;

  /** \brief Linear solver
    * Create a solver for linear systems of equations
    * Solves the linear system A*X = B or A^T*X = B for X
    * with A square and non-singular
    *
    *  If A is structurally singular, an error will be thrown during init.
    *  If A is numerically singular, the prepare step will fail.

      \generalsection{Linsol}
      \pluginssection{Linsol}

      \author Joel Andersson
      \date 2011-2016
  */
  class CASADI_EXPORT Linsol
    : public SharedObject,
      public SWIG_IF_ELSE(PrintableCommon, Printable<Linsol>) {
  public:
    /** \brief Get type name */
    static std::string type_name() {return "Linsol";}

    /// Default constructor
    Linsol();

    /// Constructor
    explicit Linsol(const std::string& name, const std::string& solver,
                    const Sparsity& sp, const Dict& opts=Dict());

    /// Access functions of the node
    LinsolInternal* operator->();
    const LinsolInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool test_cast(const SharedObjectInternal* ptr);

    /// Check if a plugin is available
    static bool has_plugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void load_plugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Query plugin name
    std::string plugin_name() const;

    /// Solve numerically
    DM solve(const DM& A, const DM& B, bool tr=false) const;

    /// Create a solve node
    MX solve(const MX& A, const MX& B, bool tr=false) const;

#ifndef SWIG
    // Symbolic factorization of the linear system, e.g. selecting pivots
    int sfact(const double* A, int mem=0) const;

    // Numeric factorization of the linear system
    int nfact(const double* A, int mem=0) const;

    // Solve factorized linear system of equations
    int solve(const double* A, double* x, int nrhs=1, bool tr=false, int mem=0) const;

#endif // SWIG
    /** \brief Number of negative eigenvalues
      * Not available for all solvers
      */
    int neig() const;

    /** \brief Matrix rank
      * Not available for all solvers
      */
    int rank() const;
  };

  /// Check if a particular plugin is available
  CASADI_EXPORT bool has_linsol(const std::string& name);

  /// Explicitly load a plugin dynamically
  CASADI_EXPORT void load_linsol(const std::string& name);

  /// Get the documentation string for a plugin
  CASADI_EXPORT std::string doc_linsol(const std::string& name);

} // namespace casadi

#endif // CASADI_LINSOL_HPP
