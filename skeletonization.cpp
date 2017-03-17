#include "skeletonization.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;
typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             Skeleton;
typedef Skeleton::vertex_descriptor                           Skeleton_vertex;
typedef Skeleton::edge_descriptor                             Skeleton_edge;

namespace {

struct Display_polylines{
  const Skeleton& skeleton;
  std::ofstream& out;
  int polyline_size;
  std::stringstream sstr;
  Display_polylines(const Skeleton& skeleton, std::ofstream& out)
    : skeleton(skeleton), out(out)
  {}
  void start_new_polyline(){
    polyline_size=0;
    sstr.str("");
    sstr.clear();
  }

  void add_node(Skeleton_vertex v){
    ++polyline_size;
    sstr << " " << skeleton[v].point;
  }

  void end_polyline()
  {
    out << polyline_size << sstr.str() << "\n";
  }
};

}

namespace skeleton {

void skeletonization(const std::string & in_file, const std::string & out_file) {
	std::ifstream input(in_file);
	if (!input) std::cout << "Failed opening in_file" << std::endl;
	Polyhedron tmesh;
	input >> tmesh;

	// extract the skeleton
	Skeleton skeleton;
	CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
	std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
	std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";

	std::ofstream output(out_file);
	BOOST_FOREACH(Skeleton_vertex v, vertices(skeleton)) {
		output << v + 1 << " " << skeleton[v].point << std::endl;
	}
	output << "#" << std::endl;
	BOOST_FOREACH(Skeleton_edge e, edges(skeleton)) {
		auto v = source(e, skeleton);
		auto u = target(e, skeleton);
		output << v + 1 << " " << u + 1 << std::endl;
	}
	output << "#" << std::endl;
	output.close();

	//Display_polylines display(skeleton,output);
	//CGAL::split_graph_into_polylines(skeleton, display);
}

} // namespace skeleton