// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "arap_support.h"
#include "igl/colon.h"
#include "igl/cotmatrix.h"
#include "igl/massmatrix.h"
#include "igl/group_sum_matrix.h"
#include "igl/covariance_scatter_matrix.h"
#include "igl/speye.h"
#include "igl/mode.h"
#include "igl/project_isometrically_to_plane.h"
#include "igl/slice.h"
#include "igl/arap_rhs.h"
#include "igl/repdiag.h"
#include "igl/columnize.h"
#include "igl/polar_svd3x3.h"
#include "igl/repmat.h"
#include "igl/verbose.h"
#include "igl/polar_dec.h"
#include "igl/polar_svd.h"
#include "igl/C_STR.h"
#include "igl/triangle_triangle_adjacency.h"

#include <cassert>
#include <iostream>

namespace igl {

	template <typename DerivedS, typename DerivedD>
	IGL_INLINE void
		fit_rotations(
			const Eigen::PlainObjectBase<DerivedS> & S,
			const bool single_precision,
			const Eigen::PlainObjectBase<DerivedD> & oldR,
			const Eigen::PlainObjectBase<DerivedD> & FArea,
			ARAPSupportData & data,
			Eigen::PlainObjectBase<DerivedD> & R)
	{
		using namespace std;
		const int dim = S.cols();
		const int nr = S.rows() / dim;
		assert(nr * dim == S.rows());
		assert(dim == 3);

		// resize output
		R.resize(dim, dim*nr); // hopefully no op (should be already allocated)

							   //std::cout<<"S=["<<std::endl<<S<<std::endl<<"];"<<std::endl;
							   //MatrixXd si(dim,dim);
		Eigen::Matrix<typename DerivedS::Scalar, 3, 3> si;// = Eigen::Matrix3d::Identity();
														  // loop over number of rotations we're computing
		//std::ofstream RVS("RVS.txt");
		for (int r = 0; r < nr; r++)
		{
			// build this covariance matrix
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					si(i, j) = S(i*nr + r, j);
				}
			}

			si = 2 * si;
			// Smooth Rotation Enhanced As-Rigid-As-Possible Mesh Animation, TVCG2015
			Eigen::Matrix<typename DerivedS::Scalar, 3, 3> sri= 4 * 0.01* FArea(r,0) * oldR.block(0, data.face_adjacency(r, 0) * dim, dim, dim).transpose();
			for (unsigned fj = 1; fj < dim; ++fj) {
				//RVS << data.face_adjacency(r, fj) << " ";
				sri += 4 * 0.01*oldR.block(0, data.face_adjacency(r, fj) * dim, dim, dim).transpose();
			}
			//RVS << std::endl;

			for (unsigned fj = 0; fj < dim; ++fj) {
				const int current_index = data.face_adjacency(r, fj);
				//std::cout << oldR.block(0, current_index*dim, dim, dim).transpose() << std::endl;
				si = si+ 4 * 0.01*oldR.block(0, current_index*dim, dim, dim).transpose();
			}
			
			typedef Eigen::Matrix<typename DerivedD::Scalar, 3, 3> Mat3;
			typedef Eigen::Matrix<typename DerivedD::Scalar, 3, 1> Vec3;
			Mat3 ri;
			if (single_precision)
			{
				polar_svd3x3(si, ri);
			}
			else
			{
				Mat3 ti, ui, vi;
				Vec3 _;
				//RVS << sri << std::endl;
				//RVS << si << std::endl;
				polar_svd(si, ri, ti, ui, _, vi);
				//RVS << ri << std::endl << std::endl;
			}
			assert(ri.determinant() >= 0);
			R.block(0, r*dim, dim, dim) = ri.block(0, 0, dim, dim).transpose();
			//cout<<matlab_format(si,C_STR("si_"<<r))<<endl;
			//cout<<matlab_format(ri.transpose().eval(),C_STR("ri_"<<r))<<endl;
		}
		//RVS.close();
	}

	template <
		typename DerivedV,
		typename DerivedF,
		typename Derivedb>
		IGL_INLINE bool arap_precomputation(
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			const int dim,
			const Eigen::PlainObjectBase<Derivedb> & b,
			ARAPSupportData & data)
	{
		using namespace std;
		using namespace Eigen;
		typedef typename DerivedV::Scalar Scalar;
		// number of vertices
		const int n = V.rows();
		data.n = n;
		assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
		assert((b.size() == 0 || b.minCoeff() >= 0) && "b out of bounds");
		// remember b
		data.b = b;
		//assert(F.cols() == 3 && "For now only triangles");
		// dimension
		//const int dim = V.cols();
		assert((dim == 3 || dim == 2) && "dim should be 2 or 3");
		data.dim = dim;
		//assert(dim == 3 && "Only 3d supported");
		// Defaults
		data.f_ext = MatrixXd::Zero(n, data.dim);

		assert(data.dim <= V.cols() && "solve dim should be <= embedding");
		bool flat = (V.cols() - data.dim) == 1;

		DerivedV plane_V;
		DerivedF plane_F;
		typedef SparseMatrix<Scalar> SparseMatrixS;
		SparseMatrixS ref_map, ref_map_dim;
		if (flat)
		{
			project_isometrically_to_plane(V, F, plane_V, plane_F, ref_map);
			repdiag(ref_map, dim, ref_map_dim);
		}
		const PlainObjectBase<DerivedV>& ref_V = (flat ? plane_V : V);
		const PlainObjectBase<DerivedF>& ref_F = (flat ? plane_F : F);
		SparseMatrixS L;
		cotmatrix(V, F, L);

		ARAPEnergyType eff_energy = data.energy;
		if (eff_energy == ARAP_ENERGY_TYPE_DEFAULT)
		{
			switch (F.cols())
			{
			case 3:
				if (data.dim == 3)
				{
					eff_energy = ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
				}
				else
				{
					eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;
				}
				break;
			case 4:
				eff_energy = ARAP_ENERGY_TYPE_ELEMENTS;
				break;
			default:
				assert(false);
			}
		}


		// Get covariance scatter matrix, when applied collects the covariance
		// matrices used to fit rotations to during optimization
		covariance_scatter_matrix(ref_V, ref_F, eff_energy, data.CSM);

		triangle_triangle_adjacency(F, data.face_adjacency);

		if (flat)
		{
			data.CSM = (data.CSM * ref_map_dim.transpose()).eval();
		}
		assert(data.CSM.cols() == V.rows()*data.dim);

		// Get group sum scatter matrix, when applied sums all entries of the same
		// group according to G
		SparseMatrix<double> G_sum;
		if (data.G.size() == 0)
		{
			if (eff_energy == ARAP_ENERGY_TYPE_ELEMENTS)
			{
				speye(F.rows(), G_sum);
			}
			else
			{
				speye(n, G_sum);
			}
		}
		else
		{
			// groups are defined per vertex, convert to per face using mode
			if (eff_energy == ARAP_ENERGY_TYPE_ELEMENTS)
			{
				Eigen::Matrix<int, Eigen::Dynamic, 1> GG;
				MatrixXi GF(F.rows(), F.cols());
				for (int j = 0; j < F.cols(); j++)
				{
					Matrix<int, Eigen::Dynamic, 1> GFj;
					slice(data.G, F.col(j), GFj);
					GF.col(j) = GFj;
				}
				mode<int>(GF, 2, GG);
				data.G = GG;
			}
			//printf("group_sum_matrix()\n");
			group_sum_matrix(data.G, G_sum);
		}
		SparseMatrix<double> G_sum_dim;
		repdiag(G_sum, data.dim, G_sum_dim);	// 沿着对角线重复G_sum data.dim次 -> 3#F x 3#F
		assert(G_sum_dim.cols() == data.CSM.rows());
		data.CSM = (G_sum_dim * data.CSM).eval();

		arap_rhs(ref_V, ref_F, data.dim, eff_energy, data.K);
		if (flat)
		{
			data.K = (ref_map_dim * data.K).eval();
		}
		assert(data.K.rows() == data.n*data.dim);

		SparseMatrix<double> Q = (-L).eval();

		if (data.with_dynamics)
		{
			const double h = data.h;
			assert(h != 0);
			SparseMatrix<double> M;
			massmatrix(V, F, MASSMATRIX_TYPE_DEFAULT, data.M);
			const double dw = (1. / data.ym)*(h*h);
			SparseMatrix<double> DQ = dw * 1. / (h*h)*data.M;
			Q += DQ;
			// Dummy external forces
			data.f_ext = MatrixXd::Zero(n, data.dim);
			data.vel = MatrixXd::Zero(n, data.dim);
		}
		// data.solver_data.solver_type = igl::min_quad_with_fixed_data<double>::QR_LLT;
		return min_quad_with_fixed_precompute(
			Q, b, SparseMatrix<double>(), true, data.solver_data);
	}

	template <
		typename Derivedbc,
		typename DerivedU,
		typename DerivedD>
		IGL_INLINE bool arap_solve(
			const Eigen::PlainObjectBase<Derivedbc> & bc,
			const Eigen::PlainObjectBase<DerivedD> & FArea,
			ARAPSupportData & data,
			Eigen::PlainObjectBase<DerivedU> & U)
	{
		using namespace Eigen;
		using namespace std;
		assert(data.b.size() == bc.rows());
		if (bc.size() > 0)
		{
			assert(bc.cols() == data.dim && "bc.cols() match data.dim");
		}
		const int n = data.n;
		int iter = 0;
		if (U.size() == 0)
		{
			// terrible initial guess.. should at least copy input mesh
#ifndef NDEBUG
			cerr << "arap_solve: Using terrible initial guess for U. Try U = V." << endl;
#endif
			U = MatrixXd::Zero(data.n, data.dim);
		}
		else
		{
			assert(U.cols() == data.dim && "U.cols() match data.dim");
		}
		// changes each arap iteration
		MatrixXd U_prev = U;
		// doesn't change for fixed with_dynamics timestep
		MatrixXd U0;
		if (data.with_dynamics)
		{
			U0 = U_prev;
		}

		Eigen::MatrixXd r_identity = Eigen::MatrixXd::Identity(3,3);
		MatrixXd oldR;
		repmat(r_identity, 1, FArea.rows(), oldR);

		while (iter < data.max_iter)
		{
			U_prev = U;
			// enforce boundary conditions exactly
			for (int bi = 0; bi < bc.rows(); bi++)
			{
				U.row(data.b(bi)) = bc.row(bi);
			}

			const auto & Udim = U.replicate(data.dim, 1);

			//std::cout << "ONE: " << U.rows() << " " << U.cols() << std::endl;
			//std::cout << "TWO: " << Udim.rows() << " " << Udim.cols() << std::endl;

			assert(U.cols() == data.dim);
			// As if U.col(2) was 0
			MatrixXd S = data.CSM * Udim;
			// THIS NORMALIZATION IS IMPORTANT TO GET SINGLE PRECISION SVD CODE TO WORK
			// CORRECTLY.
			S /= S.array().abs().maxCoeff();

			const int Rdim = data.dim;
			MatrixXd R(Rdim, data.CSM.rows());
			if (R.rows() == 2)
			{
				fit_rotations(S, true, oldR, FArea, data, R);
			}
			else
			{
				// R = 3 x 3#V or 3#F
				fit_rotations(S, false, oldR, FArea, data, R);
				//std::cout << "S: " << data.CSM.rows() << " " << data.CSM.cols() << std::endl;
				//#ifdef __SSE__ // fit_rotations_SSE will convert to float if necessary
				//      fit_rotations_SSE(S,R);
				//#else
				//      fit_rotations(S,true,R);
				//#endif
			}
			//std::cout << R.cols() << " x " << R.rows() << std::endl;
			//system("pause");
			//for(int k = 0;k<(data.CSM.rows()/dim);k++)
			//{
			//  R.block(0,dim*k,dim,dim) = MatrixXd::Identity(dim,dim);
			//}

			// record old rotations
			oldR = R;

			// Number of rotations: #vertices or #elements
			int num_rots = data.K.cols() / Rdim / Rdim;
			// distribute group rotations to vertices in each group
			MatrixXd eff_R;
			if (data.G.size() == 0)
			{
				// copy...
				eff_R = R;
			}
			else
			{
				eff_R.resize(Rdim, num_rots*Rdim);
				for (int r = 0; r < num_rots; r++)
				{
					eff_R.block(0, Rdim*r, Rdim, Rdim) =
						R.block(0, Rdim*data.G(r), Rdim, Rdim);
				}
			}
			VectorXd Rcol;
			columnize(eff_R, num_rots, 2, Rcol);
			VectorXd Bcol = -data.K * Rcol;
			assert(Bcol.size() == data.n*data.dim);
			for (int c = 0; c < data.dim; c++)
			{
				VectorXd Uc, Bc, bcc, Beq;
				Bc = Bcol.block(c*n, 0, n, 1);
				if (bc.size() > 0)
				{
					bcc = bc.col(c);
				}
				bool is_solved = min_quad_with_fixed_solve(
					data.solver_data,
					Bc, bcc, Beq,
					Uc);
				if (!is_solved) {
					system("pause");
				}
				U.col(c) = Uc;
			}

			iter++;
		}
		if (data.with_dynamics)
		{
			// Keep track of velocity for next time
			data.vel = (U - U0) / data.h;
		}

		return true;
	}



}
