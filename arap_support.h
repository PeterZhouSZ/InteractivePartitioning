// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_ARAP_SUPPORT_H
#define IGL_ARAP_SUPPORT_H
#include "igl/igl_inline.h"
#include "igl/min_quad_with_fixed.h"
#include "igl/ARAPEnergyType.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace igl
{
	struct ARAPSupportData
	{
		// n  #V
		// G  #V list of group indices (1 to k) for each vertex, such that vertex i
		//    is assigned to group G(i)
		// energy  type of energy to use
		// with_dynamics  whether using dynamics (need to call arap_precomputation
		//   after changing)
		// f_ext  #V by dim list of external forces
		// vel  #V by dim list of velocities
		// h  dynamics time step
		// ym  ~Young's modulus smaller is softer, larger is more rigid/stiff
		// max_iter  maximum inner iterations
		// K  rhs pre-multiplier
		// M  mass matrix
		// solver_data  quadratic solver data
		// b  list of boundary indices into V
		// dim  dimension being used for solving
		int n;
		Eigen::VectorXi G;
		ARAPEnergyType energy;
		bool with_dynamics;
		Eigen::MatrixXd f_ext, vel, w_kl;
		double h;
		double ym;
		int max_iter;
		Eigen::SparseMatrix<double> K, M;
		Eigen::SparseMatrix<double> CSM;
		min_quad_with_fixed_data<double> solver_data;
		Eigen::VectorXi b;
		int dim;
		Eigen::MatrixXi face_adjacency;
		ARAPSupportData() :
			n(0),
			G(),
			energy(ARAP_ENERGY_TYPE_DEFAULT),
			with_dynamics(false),
			f_ext(),
			h(1),
			ym(1),
			max_iter(10),
			K(),
			CSM(),
			solver_data(),
			b(),
			dim(-1) // force this to be set by _precomputation
		{
		};
	};

	// Compute necessary information to start using an ARAP deformation
	//
	// Inputs:
	//   V  #V by dim list of mesh positions
	//   F  #F by simplex-size list of triangle|tet indices into V
	//   dim  dimension being used at solve time. For deformation usually dim =
	//     V.cols(), for surface parameterization V.cols() = 3 and dim = 2
	//   b  #b list of "boundary" fixed vertex indices into V
	// Outputs:
	//   data  struct containing necessary precomputation
	template <
		typename DerivedV,
		typename DerivedF,
		typename Derivedb>
		IGL_INLINE bool arap_precomputation(
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			const int dim,
			const Eigen::PlainObjectBase<Derivedb> & b,
			ARAPSupportData & data);
	// Inputs:
	//   bc  #b by dim list of boundary conditions
	//   data  struct containing necessary precomputation and parameters
	//   U  #V by dim initial guess
	template <
		typename Derivedbc,
		typename DerivedU,
		typename DerivedD>
		IGL_INLINE bool arap_solve(
			const Eigen::PlainObjectBase<Derivedbc> & bc,
			const Eigen::PlainObjectBase<DerivedD> & FArea,
			ARAPSupportData & data,
			Eigen::PlainObjectBase<DerivedU> & U);

};

template bool igl::arap_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, igl::ARAPSupportData&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&);
template bool igl::arap_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, igl::ARAPSupportData&);
#endif
