#include <sstream>
#include <twg/image/image_drawing.h>

#include "logic.h"
#include "objects.h"
#include "find_borders.h"

using namespace twg;

//-----------------------------------------------------------------------------
Polygon_d calc_circle(circle a, const FindBorders& brd) {
	vec2 pos = brd.toImg(vec2(a.c.x, a.c.y));
	double r = brd.toImg(a.r);
	return computeEllipse({r, r}).move(Point_d(pos.x, pos.y));
}

//-----------------------------------------------------------------------------
class Picture
{
public:
	virtual void init_brd(FindBorders& brd) const {
		brd.process(vec2(0, 0));
		for (int i = 0; i < m_solved.process.size(); ++i)
			brd.process(vec2(m_solved.process[i][0], m_solved.process[i][1]));
	}
	virtual double getResidual(vec2 a) const = 0;
	virtual void draw(ImageDrawing_aa& img, const FindBorders& brd) const = 0;
	vector<vec2> getProcess(void) const {
		vector<vec2> result;
		for (auto& i : m_solved.process)
			result.push_back(vec2(i[0], i[1]));
		return result;
	}
protected:
	solved_t m_solved;
	fnm_f m_f;
};

//-----------------------------------------------------------------------------
class Picture_three_circles : public Picture
{
public:
	Picture_three_circles(circle a, circle b, circle c, xn_t start) : m_a(a), m_b(b), m_c(c), m_res_circle(a) {
		auto p = three_circles(m_a, m_b, m_c, false, true, true);
		m_f = p.second;
		m_solved = solve(p.first, p.second, start, 100, 1e-10, false);
		m_res_circle.c = point(m_solved.point[0], m_solved.point[1]);
		m_res_circle.r = m_solved.point[2];
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		brd.process(m_c);
		brd.process(m_res_circle);
		Picture::init_brd(brd);
	}

	double getResidual(vec2 a) const {
		double r = fabs(sqrt(sqr(m_a.c.x-a.x)+sqr(m_a.c.y-a.y))-m_a.r);
		return length(calc_vector_function(m_f, {a.x, a.y, r}));
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));
		img.drawPolyline(calc_circle(m_a, brd));
		img.drawPolyline(calc_circle(m_b, brd));
		img.drawPolyline(calc_circle(m_c, brd));

		img.setPen(Pen(1, Orange));
		img.drawPolyline(calc_circle(m_res_circle, brd));
	}
private:
	circle m_a, m_b, m_c, m_res_circle;
};

//-----------------------------------------------------------------------------
class Picture_three_lines : public Picture
{
public:
	Picture_three_lines(vec2 a, vec2 b, vec2 c, xn_t start) : m_a(a), m_b(b), m_c(c) {
		point pa(a.x, a.y), pb(b.x, b.y), pc(c.x, c.y);
		auto p = three_lines({pa, pb}, {pb, pc}, {pa, pc});
		m_f = p.second;
		m_solved = solve(square_cast_2(get_sle_function(calc_jacobi_matrix_numeric_functon(p.second), p.second)), p.second, start, 100, 1e-10, false);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		brd.process(m_c);
		Picture::init_brd(brd);
	}

	double getResidual(vec2 a) const {
		return length(calc_vector_function(m_f, {a.x, a.y}));
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));

		auto draw_line = [&] (vec2 a, vec2 b) {
			vec2 start, end;
			start = a + (b-a) * (-5);
			end =   a + (b-a) * (5);
			img.drawLine(brd.toImg(start), brd.toImg(end));
		};

		draw_line(m_a, m_b);
		draw_line(m_a, m_c);
		draw_line(m_b, m_c);
	}
private:
	vec2 m_a, m_b, m_c;
};

//-----------------------------------------------------------------------------
class Picture_sin_and_line : public Picture
{
public:
	Picture_sin_and_line(vec2 a, xn_t start) : m_a(a), m_b(0, 0) {
		point pa(a.x, a.y);
		auto p = sin_and_line(pa);
		m_f = p.second;
		m_solved = solve(p.first, p.second, start, 100, 1e-10, false);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		brd.process(vec2(-3, -3));
		Picture::init_brd(brd);
	}

	double getResidual(vec2 a) const {
		return length(calc_vector_function(m_f, {a.x, a.y}));
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));

		auto draw_line = [&] (vec2 a, vec2 b) {
			vec2 start, end;
			start = a + (b-a) * (-5);
			end =   a + (b-a) * (5);
			img.drawLine(brd.toImg(start), brd.toImg(end));
		};

		draw_line(m_a, m_b);

		vec2 last;
		for (int i = 0; i < 1000; ++i) {
			double x = (i - 500.0) / 100.0;
			double y = sin(x);
			vec2 current(x, y);
			if (i != 0)
				img.drawLine(brd.toImg(last), brd.toImg(current));
			last = current;
		}
	}
private:
	vec2 m_a, m_b;
};


//-----------------------------------------------------------------------------
void draw_picture(const Picture& pic, wstring filename) {
	FindBorders brd(500, 20, true);
	pic.init_brd(brd);
	brd.finish();

	// Создание изображения
	Point_i size(brd.getCalculatedSize().x, brd.getCalculatedSize().y);
	ImageDrawing_aa img(size);

	vector<vector<double>> mas(size.x, vector<double>(size.y));
	img.clear(White);

	// Заполняем массив невязкой по всем координатам
	double max1 = -100000000000;
	double min1 = 100000000000;
	for (int i = 0; i < size.x; i++) {
		for (int j = 0; j < size.y; j++) {
			mas[i][j] = pic.getResidual(brd.fromImg(vec2(i, j)));
			min1 = min(min1, mas[i][j]);
			max1 = max(max1, mas[i][j]);
		}
	}

	// Рисуем все нормированные точки из массива
	for (int i = 0; i < size.x; i++) {
		for (int j = 0; j < size.y; j++) {
			double v = mas[i][j];
			double pos = (v - min1)/(max1-min1);
			pos = int(pos*20)/20.0;
			img[Point_i(i, j)] = getColorBetween(pos, White, Black);
		}
	}

	// Рисуем оси координат
	img.setPen(Pen(1.5, Red));
	vec2 xa = brd.toImg(vec2(-1, 0));
	vec2 xb = brd.toImg(vec2(1, 0));
	vec2 ya = brd.toImg(vec2(0, -1));
	vec2 yb = brd.toImg(vec2(0, 1));
	if (xa.y >= 0 && xa.y <= size.y) {
		xa.x = -1;
		xb.x = size.x;
		img.drawLine(xa, xb);
	} 
	if (ya.x >= 0 && ya.x <= size.x) {
		ya.y = -1;
		yb.y = size.y;
		img.drawLine(ya, yb);
	}

	pic.draw(img, brd);

	// Рисуем линии поиска решения
	auto process = pic.getProcess();
	for (int i = 0; i < process.size()-1; i++) {
		auto a = brd.toImg(process[i]);
		auto b = brd.toImg(process[i+1]);

		img.setPen(Pen(1, Green));
		img.drawLine(a, b);
	}

	// Рисуем точки поиска решения
	for (int i = 0; i < process.size(); i++) {
		auto a = brd.toImg(process[i]);

		img.setBrush(Brush(Green));
		img.drawPolygon(computeEllipse(Point_d(1.5, 1.5)).move(a));
	}

	saveToBmp(&img, filename);
}

int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
	// Инициализация рисунка
	/*circle a(0, 0, 1);
	circle b(5, 0, 1);
	circle c(10, 5, 5);
	Picture_three_circles pic(a, b, c, {-1, -1, 13});
	draw_picture(pic, L"three_circles_4lab.bmp");*/

	vec2 a(0.5, 0.5), b(3, 1), c(-1, 3);
	//vec2 a(0, 0), b(1, 2), c(2, 1); 
	Picture_three_lines pic(a, b, c, {-3, -3});
	draw_picture(pic, L"three_lines_4lab.bmp");

	/*vec2 a(0.5, 3);
	Picture_sin_and_line pic(a, {3, 5});
	draw_picture(pic, L"sin_and_line_4lab.bmp");*/
}