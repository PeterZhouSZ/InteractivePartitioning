#include "arap.h"
#include <igl/per_face_normals.h>
void deform_selected_faces(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
	Eigen::MatrixXd N;
	igl::per_face_normals(V, F, N);
	Eigen::Vector3d printDir(0, 1, 0);
	Eigen::VectorXd NRisky;

	NRisky = N * printDir;

	for (unsigned i = 0; i < NRisky.rows(); ++i) {
		if (NRisky(i) + 0.707 < 0) {

		}
	}
}