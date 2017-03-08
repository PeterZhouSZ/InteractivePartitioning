#include "interactiveSelection.h"
#include <igl/adjacency_list.h>
#include <igl/triangle_triangle_adjacency.h>
#include <iostream>
#include <queue>

#define INF 1e30
namespace intRobo {
	typedef std::pair<int, double> TupleFD;
	struct opt {
		bool operator()(const TupleFD& r1, const TupleFD& r2) {
			return r1.second > r2.second;
		}
	};

	IGL_INLINE double vert_dist(const Eigen::MatrixXd & V, const int source, const int target) {
		return (V.row(source) - V.row(target)).squaredNorm();

	}

	void part_selection(const int f, const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, std::vector<int>& selection) {
		std::vector<std::vector<int>> adjacencyList;
		igl::adjacency_list(F, adjacencyList, true);
		
		std::vector<double> vec_dist(F.rows(), INF);
		std::priority_queue<TupleFD, std::vector<TupleFD>, opt> pq;

		vec_dist[f] = 0;
		pq.push(TupleFD(f, vec_dist[f]));
		while (!pq.empty()) {
			TupleFD v = pq.top();
			pq.pop();
			if (selection.size() == 500) return;
			selection.push_back(v.first);
			for (int ng : adjacencyList[v.first]) {
				double alt = vert_dist(V, v.first, ng);
				if (vec_dist[ng] > v.second + alt) {
					vec_dist[ng] = v.second + alt;
					pq.push(TupleFD(ng, vec_dist[ng]));
				}
			}
		}
	}
}