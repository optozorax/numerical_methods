#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <twg/image/image_drawing.h>

#include "logic.h"
#include "objects.h"
#include "find_borders.h"

using namespace twg;

#define TWO_CIRCLES
//#define THREE_CIRCLES
//#define CUBIC_FUNCTION

int main() {
	FindBorders brd(500, 10, true);
	brd.process({-2, -2});
	brd.process({2, 2});
	brd.finish();

	#ifdef TWO_CIRCLES
	circle a(-0.5, 0, 1), b(0.5, 0, 1);
	auto snu = two_circles(a, b).first;
	#endif

	#ifdef THREE_CIRCLES
	circle a(-0.5, 0, 1), b(0.5, 0, 1), c(0, 0, 0.5); bool a_in = true, b_in = true, c_in = false;
	auto snu = three_circles(a, b, c, a_in, b_in, c_in).first;
	#endif

	#ifdef CUBIC_FUNCTION
	fn_f f1 = [] (const xn_t& x) -> double { return x[0]*x[0]*x[0]+x[1]*x[1]*x[1]-1; };
	fn_f f2 = [] (const xn_t& x) -> double { return (x[0]-1)*(x[0]-1)*(x[0]-1)+x[1]*x[1]-0.5; };
	fnmv_f f = {f1, f2};

	jnm_f j = [] (const xn_t& x) -> matrix_t {
		matrix_t result(2, xn_t(2));
		result[0][0] = 3.0*x[0]*x[0];
		result[0][1] = 3.0*x[1]*x[1];

		result[1][0] = 3.0*(x[0]-1)*(x[0]-1);
		result[1][1] = 2.0*x[1];

		return result;
	};

	auto snu = get_sle_function(j, f);
	#endif

	ImageDrawing_aa img(brd.getCalculatedSize());
	img.clear(White);

	vector<Color> colors = {Red, Blue, Green, Orange, Miku};
	vector<xn_t> answers;
	for (int i = 0; i < img.width(); ++i) {
		cout << "\r" << fixed << setprecision(2) << setw(6) << double(100.0 * i)/img.width() << "%";
		for (int j = 0; j < img.height(); ++j) {
			vec2 p(i, j);
			p = brd.fromImg(p);

			#ifdef TWO_CIRCLES
			xn_t x0(2);
			x0[0] = p.x; 
			x0[1] = p.y;
			#endif

			#ifdef THREEE_CIRCLES
			xn_t x0(3); 
			x0[0] = p.x; 
			x0[1] = p.y;
			x0[2] = 2.22;
			#endif 

			#ifdef CUBIC_FUNCTION
			xn_t x0(2);
			x0[0] = p.x; 
			x0[1] = p.y;
			#endif

			auto result = solve(snu, x0, 1000, 1e-9, false);
			result.point[2] = fabs(result.point[2]);
			
			if (result.residual < 1e-5) {
				int pos = -1;
				bool isFinded = false;
				for (int k = 0; k < answers.size(); k++) {
					if (length(answers[k] - result.point) < 0.01) {
						isFinded = true;
						pos = k;
						break;
					}
				}
				if (!isFinded) {
					answers.push_back(result.point);
					pos = answers.size() - 1;
					isFinded = true;
					cout << "added: " << result.point << endl;
				}
				if (pos > colors.size()) pos = colors.size();
				img[Point_i(i, j)] = colors[pos];
			}
		}
	}

	img.setPen(Pen(2, Black));
	auto draw_circle = [&img, &brd] (circle a) {
		img.drawPolyline(computeEllipse(brd.toImg({a.r, a.r}) - brd.toImg({0, 0})).move(brd.toImg({a.c.x, a.c.y})));;
	};

	#ifdef TWO_CIRCLES
	draw_circle(a);
	draw_circle(b);
	#endif

	#ifdef THREE_CIRCLES
	draw_circle(a);
	draw_circle(b);
	draw_circle(c);
	#endif

	#ifdef CUBIC_FUNCTION
	#endif

	int counter = 0;
	for (auto& i : answers) {
		img.setPen(Pen(2, getColorBetween(0.5, colors[counter], White)));
		if (i.size() == 3)
			draw_circle(circle(i[0], i[1], i[2]));
		else
			draw_circle(circle(i[0], i[1], (brd.fromImg({2, 2})-brd.fromImg({0, 0})).x));
		counter++;
	}

	saveToPng(&img, L"newton_fractal.png");
}