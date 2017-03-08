#include <igl/avg_edge_length.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/read_triangle_mesh.h>
#include <igl/embree/ambient_occlusion.h>
#include <igl/viewer/Viewer.h>
#include <igl/viewer/ViewerPlugin.h>
#include <igl/jet.h>
#include <igl/barycenter.h>

#include <nanogui/formhelper.h>
#include <nanogui/screen.h>
#include <nanogui/messagedialog.h>
#include <nanogui/slider.h>
#include <nanogui/progressbar.h>
#include <nanogui/toolbutton.h>

#include <iostream>
#include "qJet.h"

#include <tinyxml2.h>

#include <GLFW/glfw3.h>

// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::VectorXd AO;
Eigen::Vector3d vecZDir;

namespace meshProcessing {
	Eigen::Vector3d pmin, pmax;
	double coff;
	void normalizeMesh(Eigen::MatrixXi& F, Eigen::MatrixXd& V) {
		pmin = V.colwise().minCoeff();
		pmax = V.colwise().maxCoeff();

		coff = std::fabs((pmax - pmin).maxCoeff());

		for (auto i = 0; i < V.rows(); ++i) {
			V.row(i) /= coff;
		}

		// Compute the center point of the mesh
		Eigen::Vector3d center(0, 0, 0);
		Eigen::MatrixXd centers;
		igl::barycenter(V, F, centers);
		center = centers.colwise().mean();		// Mean value of barycenter matrix
		for (auto i = 0; i < V.rows(); ++i) {
			V.row(i) -= center;
		}

	}
}

// It allows to change the degree of the field when a number is pressed
bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	using namespace Eigen;
	using namespace std;
	const RowVector3d color(0.9, 0.85, 0.9);
	switch (key)
	{
	case '1':
		// Show the mesh without the ambient occlusion factor
		viewer.data.set_colors(color);
		break;
	case '2':
	{
		// Show the mesh with the ambient occlusion factor
		MatrixXd C = color.replicate(V.rows(), 1);
		for (unsigned i = 0; i < C.rows(); ++i)
			C.row(i) *= AO(i);//std::min<double>(AO(i)+0.2,1);
		viewer.data.set_colors(C);
		break;
	}
	case '.':
		viewer.core.lighting_factor += 0.1f;
		break;
	case ',':
		viewer.core.lighting_factor -= 0.1f;
		break;
	default: break;
	}
	viewer.core.lighting_factor =
		std::min(std::max(viewer.core.lighting_factor, 0.f), 1.f);

	return false;
}


void someOpt(igl::viewer::Viewer &viewer) {

	// Load a mesh in OFF format
	meshProcessing::normalizeMesh(F, V);

	Eigen::MatrixXd N;
	igl::per_vertex_normals(V, F, N);

	// Compute ambient occlusion factor using embree
	igl::embree::ambient_occlusion(V, F, V, N, 500, AO);
	AO = 1.0 - AO.array();

	return;
}

int main(int argc, char *argv[])
{
	using namespace std;
	using namespace Eigen;

	cout <<
		"Press 1 to turn off Ambient Occlusion" << endl <<
		"Press 2 to turn on Ambient Occlusion" << endl <<
		"Press . to turn up lighting" << endl <<
		"Press , to turn down lighting" << endl;

	// load external setting file
	tinyxml2::XMLDocument doc;
	if (!doc.LoadFile("settings.xml") == tinyxml2::XML_SUCCESS) {
		std::cout << "Couldn't open settings.xml, will automatically exit." << std::endl;
		return -1;
	};
	

	// show mesh
	igl::viewer::Viewer viewer;

	viewer.callback_init = [&](igl::viewer::Viewer &viewer) {
		using namespace nanogui;

		//viewer.ngui->window()->setVisible(false);
		viewer.ngui->window()->setTitle("Basic panel");
		viewer.ngui->window()->setPosition(Eigen::Vector2i(10, 80));
		// Main menu window
		viewer.ngui->setFixedSize(Eigen::Vector2i(120, 20));
		nanogui::Window *window = viewer.ngui->addWindow(Eigen::Vector2i(1080, 10), "Support-free Opts");
		viewer.ngui->addGroup("Main Menu");

		Widget *panel = new Widget(window);
		viewer.ngui->addWidget("Line width", panel);
		panel->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 20));

		Slider *slider = new Slider(panel);
		slider->setValue(0.5f);
		slider->setFixedWidth(80);
		TextBox *textBox = new TextBox(panel);
		textBox->setFixedSize(Vector2i(60, 25));
		textBox->setValue("50");
		textBox->setUnits("%");
		slider->setCallback([textBox, &viewer](float value) {
			textBox->setValue(std::to_string((int)(value * 100)));
			viewer.core.line_width = 2 * value;
		});


		panel = new Widget(window);
		viewer.ngui->addWidget("Anesthetic", panel);
		panel->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 20));

		slider = new Slider(panel);
		slider->setValue(0.5f);
		slider->setFixedWidth(80);
		textBox = new TextBox(panel);
		textBox->setFixedSize(Vector2i(60, 25));
		textBox->setValue("50");
		textBox->setUnits("%");

		viewer.ngui->addButton("Run Algorithm", [&]() {
		});

		

		window = viewer.ngui->addWindow(Vector2i(10, 10), "Main menu");
		window->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 6));
		//window->setFixedWidth(1000);

		auto *pb = new PopupButton(window, "IMPORTER", ENTYPO_ICON_EXPORT);
		auto *popup = pb->popup();
		popup->setLayout(new GroupLayout());
		pb = new PopupButton(popup, "Risky face");
		auto *popup_b = pb->popup();
		popup_b->setLayout(new GridLayout(Orientation::Horizontal, 2,
			Alignment::Middle, 15, 5));
		
		const char textPars[4][2] = { "a", "b", "c", "d" };
		TextBox * textBoxPars[4];
		for (int i = 0; i < 4; ++i) {
			new Label(popup_b, std::string(textPars[i]));
			textBoxPars[i] = new TextBox(popup_b);
			textBoxPars[i]->setEditable(true);
			textBoxPars[i]->setFixedSize(Vector2i(100, 20));
			textBoxPars[i]->setValue("0.0");
			textBoxPars[i]->setFontSize(16);
			textBoxPars[i]->setFormat("[-]?[0-9]*\\.?[0-9]+");
			vecZDir[i] = 0.0;
			textBoxPars[i]->setCallback([i](const std::string& value){
				vecZDir[i] = atof(value.c_str());
				return true;
			});
		}
		textBoxPars[1]->setValue("1.0");
		auto *b = new Button(popup_b, "Detect");
		b->setCallback([&viewer]() {
			Eigen::MatrixXd N;
			Eigen::VectorXd NRisky;
			std::cout << vecZDir << std::endl;
			std::cout << vecZDir << std::endl;
			const Eigen::RowVector3d color(0.9, 0.85, 0);
			Eigen::MatrixXd C = color.replicate(viewer.data.V.rows(), 1);
			C.resize(viewer.data.V.rows(), 3);

			igl::per_face_normals(viewer.data.V, viewer.data.F, N);
			NRisky = N * vecZDir;

			for (unsigned i = 0; i < NRisky.rows(); ++i) {
				if (NRisky(i) + 0.707 < 0) {
					auto &VRisky = viewer.data.F.row(i);
					for (unsigned j = 0; j < VRisky.cols(); ++j) {
						C.row(VRisky[j]) = Eigen::RowVector3d(1.0, 0, 0);
					}
				}
			}
			viewer.data.set_colors(C);
		});

		b = new Button(popup, "Import stress");
		//b->setFlags(Button::ToggleButton);
		b->setCallback([&viewer]() {

		});

		pb = new PopupButton(window, "WEIGHT", ENTYPO_ICON_EXPORT);
		popup = pb->popup();
		// Wights options window

		GridLayout *layout =
			new GridLayout(Orientation::Horizontal, 2,
				Alignment::Middle, 15, 5);
		layout->setColAlignment(
		{ Alignment::Maximum, Alignment::Fill });
		layout->setSpacing(0, 10);
		popup->setLayout(layout);

		std::vector<IntBox<float>* > refWeights;
		auto elem = doc.FirstChildElement("Weights")->FirstChildElement("Value");
		for (int i = 0; i < 6; ++i) {
			string strLabel = string("w_").append(std::to_string(i + 1)).append(" :");
			new Label(popup, strLabel, "sans-bold");
			auto intBox = new IntBox<float>(popup);
			intBox->setEditable(true);
			intBox->setFixedSize(Vector2i(100, 20));
			intBox->setValue(atof(elem->FirstChild()->ToText()->Value()));
			elem = elem->NextSiblingElement();
			intBox->setDefaultValue("0");
			intBox->setFontSize(16);
			intBox->setFormat("[-]?[0-9]*\\.?[0-9]+");
			intBox->setSpinnable(true);
			intBox->setMinMaxValues(0, 1);
			intBox->setValueIncrement(0.02);
			intBox->value();
			refWeights.push_back(intBox);
		}



		pb = new PopupButton(window, "Risky face 1");
		
		viewer.screen->setVisible(true);
		viewer.screen->performLayout();
		//Change the title of main window
		glfwSetWindowTitle(viewer.window, doc.FirstChildElement("Window-title")->FirstChild()->Value());
		glfwSetWindowSize(viewer.window,
			atoi(doc.FirstChildElement("Window-size")->FirstChildElement("Width")->FirstChild()->Value()),
			atoi(doc.FirstChildElement("Window-size")->FirstChildElement("Height")->FirstChild()->Value()));
		return false;
	};

	//viewer.callback_key_down = &key_down;
	//key_down(viewer, '2', 0);
	viewer.core.background_color = Eigen::Vector4f(0.42f, 0.43f, 1.0f, 1.0f);
	viewer.core.show_lines = false;
	viewer.core.lighting_factor = 0.6f;
	
	viewer.launch();
}
