#include <twg/image/image_drawing.h>

#include "logic.h"
#include "objects.h"

using namespace twg;

Polygon_d calc_circle(circle a, double mx, Point_d offset) {
	return computeEllipse({a.r*mx, a.r*mx}).move(Point_d(a.c.x, a.c.y) * mx + offset);
}

int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
	circle ca(0, 0, 1);
	circle cb(5, 0, 1);
	circle cc(10, 5, 5);
	auto p = three_circles(ca, cb, cc, false, true, true);
	xn_t x_0 = {0, 0, 13};
	auto solved = solve(calc_jacobi_matrix_numeric_functon(p.first), p.first, x_0, 100, 1e-10, true);
	Point_i size(500, 500);
	Point_d offset(100, 50);
	double mx = 20;


	ImageDrawing_aa img(size);

	vector<vector<double>> mas(size.y, vector<double>(size.x));
	img.clear(White);

	// Заполняем массив невязкой по всем координатам
	double max1 = -100000000000;
	double min1 = 100000000000;
	for (int i = 0; i < size.y; i++) {
		for (int j = 0; j < size.x; j++) {
			Point_d a(i, j);
			a -= offset;
			a /= mx;
			double r = fabs(sqrt(sqr(ca.c.x-a.x)+sqr(ca.c.y-a.y))-ca.r);
			double v = length(calc(p.first, {a.x, a.y, r}));
			mas[i][j] = v;
			min1 = min(min1, mas[i][j]);
			max1 = max(max1, mas[i][j]);
		}
	}

	// Рисуем все нормированные точки из массива
	for (int i = 0; i < size.y; i++) {
		for (int j = 0; j < size.x; j++) {
			double v = mas[i][j];
			double pos = (v - min1)/(max1-min1);
			pos = int(pos*20)/20.0;
			img[Point_i(i, j)] = getColorBetween(pos, White, Black);
		}
	}

	// Рисуем линии поиска решения
	for (int i = 0; i < solved.process.size()-1; i++) {
		Point_d a(solved.process[i][0], solved.process[i][1]);
		Point_d b(solved.process[i+1][0], solved.process[i+1][1]);

		a = a*mx + offset;
		b = b*mx + offset;

		img.setPen(Pen(0.5, Green));
		img.drawLine(a, b);
	}

	// Рисуем точки поиска решения
	for (int i = 0; i < solved.process.size(); i++) {
		Point_d a(solved.process[i][0], solved.process[i][1]);
		
		a = a*mx + offset;

		img.setBrush(Brush(Green));
		img.drawPolygon(computeEllipse(Point_d(1.5, 1.5)).move(a));
	}

	img.setPen(Pen(1, Blue));
	img.drawPolyline(calc_circle(ca, mx, offset));
	img.drawPolyline(calc_circle(cb, mx, offset));
	img.drawPolyline(calc_circle(cc, mx, offset));

	img.setPen(Pen(1, Orange));
	img.drawPolyline(calc_circle({solved.point[0], solved.point[1], fabs(solved.point[2])}, mx, offset));

	saveToBmp(&img, L"4lab.bmp");
}