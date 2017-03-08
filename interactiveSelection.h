#pragma once
#include <igl/adjacency_matrix.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

namespace intRobo {
	void part_selection(
			const int f,
			const Eigen::MatrixXd & V,
			const Eigen::MatrixXi & F,
			std::vector<int>& selection);
}