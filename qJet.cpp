#include "qJet.h"
void qJetColor(double v, double color[3], double vmin, double vmax)
{
	if (v < vmin) v = vmin;
	if (v > vmax) v = vmax;
	const double dv = vmax - vmin;

	double r = 1, g = 1, b = 1;
	if (v < (vmin + 0.25 * dv)) {
		r = 0;
		g = 4 * (v - vmin) / dv;
	}
	else if (v < (vmin + 0.5 * dv)) {
		r = 0;
		b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
	}
	else if (v < (vmin + 0.75 * dv)) {
		r = 4 * (v - vmin - 0.5 * dv) / dv;
		b = 0;
	}
	else {
		g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
		b = 0;
	}
	color[0] = r;
	color[1] = g;
	color[2] = b;
	return;
}