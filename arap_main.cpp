#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/doublearea.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

#include "tutorial_shared_path.h"

#include "arap_support.h"

typedef
std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
RotationList;

const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
std::vector<int> vecBC;
Eigen::MatrixXd V, U, Area;
Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
igl::ARAPSupportData arap_data;

bool pre_draw(igl::viewer::Viewer & viewer)
{
	using namespace Eigen;
	using namespace std;
	if (!viewer.core.is_animating) {
		return false;
	}

	MatrixXd bc(b.size(), V.cols());
	for (int i = 0; i < b.size(); i++)
	{
		bc.row(i) = V.row(b(i));
		if (bc(i, 1) < 0.2)continue;
		bc(i, 0) -= 0.051;
		bc(i, 1) -= 0.051;
	}
	igl::arap_solve(bc, Area, arap_data, U);
	viewer.data.set_vertices(U);
	viewer.data.compute_normals();
	if (viewer.core.is_animating)
	{
		anim_t += anim_t_dir;
	}
	viewer.core.is_animating = false;// true; // false;
	return false;
}



bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int mods)
{
	switch (key)
	{
	case ' ':
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}

int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

#ifdef SR_ARAP
	igl::readOFF(TUTORIAL_SHARED_PATH "/bunny_smooth.off", V, F);
	U = V;
	igl::readDMAT(TUTORIAL_SHARED_PATH "/decimated-knight-selection.dmat", S);
#else
	igl::readOFF(TUTORIAL_SHARED_PATH "/decimated-knight.off", V, F);
	U = V;
	igl::readDMAT(TUTORIAL_SHARED_PATH "/decimated-knight-selection.dmat", S);
#endif

	// compute each face's area
	igl::doublearea(V, F, Area);
	std::cout << "AREA = " << Area.rows() << " " << Area.cols() << std::endl;

#ifdef SR_ARAP
	// vertices in selection
	std::vector<int> vecBoundaryFaces;
	for (unsigned i = 0; i < F.rows(); ++i) {
		if (V(F(i,0), 1) < -0.35 || V(F(i, 1), 1) < -0.35 || V(F(i, 2), 1) < -0.35) {
			vecBoundaryFaces.push_back(i);
		}
	}
	for (unsigned i = 0; i < V.rows(); ++i) {
		if (V(i, 1) < -0.35) vecBC.push_back(i);
	}

	b = Eigen::VectorXi(vecBC.size());
	for (unsigned i = 0; i < vecBC.size(); ++i) {
		b[i] = vecBC[i];
	}

	// Centroid
	mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
	// pre-computation
	// Set color based on selection
	MatrixXd C(F.rows(), 3);
	RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
	RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
	for (int f = 0; f < F.rows(); f++) {
		C.row(f) = gold;
	}

	for (int f = 0; f < vecBoundaryFaces.size(); ++f) {
		C.row(vecBoundaryFaces[f]) = purple;
	}

	arap_data.energy = igl::ARAP_ENERGY_TYPE_ELEMENTS;
#else
	// vertices in selection
	igl::colon<int>(0, V.rows() - 1, b);
	b.conservativeResize(stable_partition(b.data(), b.data() + b.size(),
		[](int i)->bool {return S(i) >= 0; }) - b.data());
	// Centroid
	mid = 0.5*(V.colwise().maxCoeff() + V.colwise().minCoeff());
	// pre-computation
	MatrixXd C(F.rows(), 3);
	RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
	RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
	for (int f = 0; f < F.rows(); f++)
	{
		if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
		{
			C.row(f) = purple;
		}
		else
		{
			C.row(f) = gold;
		}
	}
	arap_data.energy = igl::ARAP_ENERGY_TYPE_ELEMENTS;
#endif

	arap_data.max_iter = 100;
	igl::arap_precomputation(V, F, V.cols(), b, arap_data);

	// Plot the mesh with pseudo-colors
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(U, F);
	viewer.data.set_colors(C);
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.show_lines = false;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;
	cout <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}