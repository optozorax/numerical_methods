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
void draw_line(ImageDrawing_aa& img, const FindBorders& brd, vec2 a, vec2 b) {
	vec2 start, end;
	start = a + (b-a) * (-5);
	end =   a + (b-a) * (5);
	img.drawLine(brd.toImg(start), brd.toImg(end));
}

//-----------------------------------------------------------------------------
class Picture
{
public:
	void init_solve(const pair<sle_f, fnm_f>& p, const xn_t& x_0, bool is_numeric, const sqr_f& square_cast, int iterations, double residual) {
		m_f = p.second;
		if (is_numeric) {
			auto j = calc_jacobi_matrix_numeric_functon(m_f);
			auto sle = get_sle_function(j, m_f);
			m_solved = solve(square_cast(sle), p.second, x_0, iterations, residual, false);
		} else
			m_solved = solve(square_cast(p.first), p.second, x_0, iterations, residual, false);
	}
	virtual void init_brd(FindBorders& brd) const {
		brd.process(vec2(0, 0));
		for (int i = 0; i < m_solved.x_process.size(); ++i)
			brd.process(vec2(m_solved.x_process[i][0], m_solved.x_process[i][1]));
	}
	virtual double getResidual(vec2 a) const {
		return length(calc_vector_function(m_f, {a.x, a.y}));
	}
	virtual void draw(ImageDrawing_aa& img, const FindBorders& brd) const = 0;
	vector<vec2> getProcess(void) const {
		vector<vec2> result;
		for (auto& i : m_solved.x_process)
			result.push_back(vec2(i[0], i[1]));
		return result;
	}
protected:
	solved_t m_solved;
	fnm_f m_f;
};

//-----------------------------------------------------------------------------
class Picture_two_circles : public Picture
{
public:
	Picture_two_circles(circle a, circle b, xn_t x_0, bool is_numeric = false, const sqr_f& square_cast = square_cast_none, int iterations = 100, double residual = 1e-10) : m_a(a), m_b(b) {
		init_solve(two_circles(m_a, m_b), x_0, is_numeric, square_cast, iterations, residual);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		Picture::init_brd(brd);
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));
		img.drawPolyline(calc_circle(m_a, brd));
		img.drawPolyline(calc_circle(m_b, brd));
	}
private:
	circle m_a, m_b;
};

//-----------------------------------------------------------------------------
class Picture_one_circle : public Picture
{
public:
	Picture_one_circle(circle a, xn_t x_0, bool is_numeric = false, const sqr_f& square_cast = square_cast_none, int iterations = 100, double residual = 1e-10) : m_a(a) {
		init_solve(one_circle(m_a), x_0, is_numeric, square_cast, iterations, residual);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		Picture::init_brd(brd);
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));
		img.drawPolyline(calc_circle(m_a, brd));
	}
private:
	circle m_a;
};

//-----------------------------------------------------------------------------
class Picture_two_circles_and_line : public Picture
{
public:
	Picture_two_circles_and_line(circle a, circle b, line c, xn_t x_0, bool is_numeric = false, const sqr_f& square_cast = square_cast_none, int iterations = 100, double residual = 1e-10) : m_a(a), m_b(b), m_la(c.a.x, c.a.y), m_lb(c.b.x, c.b.y) {
		init_solve(two_circles_and_line(m_a, m_b, c), x_0, is_numeric, square_cast, iterations, residual);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		Picture::init_brd(brd);
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));
		img.drawPolyline(calc_circle(m_a, brd));
		img.drawPolyline(calc_circle(m_b, brd));
		draw_line(img, brd, m_la, m_lb);
	}
private:
	circle m_a, m_b;
	vec2 m_la, m_lb;
};

//-----------------------------------------------------------------------------
class Picture_three_circles : public Picture
{
public:
	Picture_three_circles(circle a, circle b, circle c, bool a_in, bool b_in, bool c_in, xn_t x_0, bool is_numeric = false, const sqr_f& square_cast = square_cast_none, int iterations = 100, double residual = 1e-10) : m_a(a), m_b(b), m_c(c), m_res_circle(a) {
		init_solve(three_circles(m_a, m_b, m_c, a_in, b_in, c_in), x_0, is_numeric, square_cast, iterations, residual);
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
	Picture_three_lines(vec2 a, vec2 b, vec2 c, xn_t x_0, bool is_numeric = false, const sqr_f& square_cast = square_cast_none, int iterations = 100, double residual = 1e-10) : m_a(a), m_b(b), m_c(c) {
		point pa(a.x, a.y), pb(b.x, b.y), pc(c.x, c.y);
		init_solve(three_lines({pa, pb}, {pb, pc}, {pa, pc}), x_0, is_numeric, square_cast, iterations, residual);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		brd.process(m_c);
		Picture::init_brd(brd);
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));		

		draw_line(img, brd, m_a, m_b);
		draw_line(img, brd, m_a, m_c);
		draw_line(img, brd, m_b, m_c);
	}
private:
	vec2 m_a, m_b, m_c;
};

//-----------------------------------------------------------------------------
class Picture_sin_and_line : public Picture
{
public:
	Picture_sin_and_line(vec2 a, xn_t x_0, bool is_numeric = false, const sqr_f& square_cast = square_cast_none, int iterations = 100, double residual = 1e-10) : m_a(a), m_b(0, 0) {
		init_solve(sin_and_line(point(a.x, a.y)), x_0, is_numeric, square_cast, iterations, residual);
	}

	void init_brd(FindBorders& brd) const {
		brd.process(m_a);
		brd.process(m_b);
		brd.process(vec2(-3, -3));
		Picture::init_brd(brd);
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.setPen(Pen(1, Blue));

		draw_line(img, brd, m_a, m_b);

		// Рисуем синусоиду
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
void draw_picture(const Picture& pic, wstring filename, int pic_size) {
	const double step_count = 20;
	FindBorders brd(pic_size - 40, 20, true);
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
			pos = int(pos*step_count)/step_count;
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
	#ifdef _DEBUG
		int pic_size = 50;
	#else
		int pic_size = 500;
	#endif

	// Инициализация рисунка
	{
	circle a(0, 0, 1), b(5, 0, 1), c(10, 5, 5);
	Picture_three_circles pic(a, b, c, true, true, true, {-1, -1, 13}, false, square_cast_none, 100, 1e-10);
	draw_picture(pic, L"three_circles_4lab.bmp", pic_size);
	}

	{
	//vec2 a(0.5, 0.5), b(3, 1), c(-1, 3);
	vec2 a(0, 0), b(1, 2), c(3, 1); 
	Picture_three_lines pic(a, b, c, {-3, -3}, false, square_cast_2, 100, 1e-10);
	draw_picture(pic, L"three_lines_4lab.bmp", pic_size);
	}

	{
	vec2 a(0.5, 3);
	Picture_sin_and_line pic(a, {3, 5}, false, square_cast_none, 100, 1e-10);
	draw_picture(pic, L"sin_and_line_4lab.bmp", pic_size);
	}

	{
	//circle a(0, 0, 1), b(1, 0.7, 1.2); // Пересекаются в двух точках
	circle a(0, 0, 0.8), b(sqrt(2.0)/2, sqrt(2.0)/2, 0.2); // Пересекаются в одной точке
	//circle a(0, 0, 0.8), b(1, 0.7, 0.3); // Не пересекаются
	Picture_two_circles pic(a, b, {-1, 0.5}, false, square_cast_none, 100, 1e-10);
	draw_picture(pic, L"two_circles_4lab.bmp", pic_size);
	}

	{
	circle a(0, 0, 0.8), b(sqrt(2.0)/2, sqrt(2.0)/2, 0.2);
	line c = {point(sqrt(2.0)/2, sqrt(2.0)/2), point(-1, 0)};
	line c = {point(sqrt(2.0)*(1.0/2.0-0.1), sqrt(2.0)*(1.0/2.0-0.1)), point(-1, 0)};
	Picture_two_circles_and_line pic(a, b, c, {-1, 0.5}, false, square_cast_3, 100, 1e-10);
	draw_picture(pic, L"two_circles_and_line_4lab.bmp", pic_size);
	}

	{
	circle a(0, 0, 1);
	Picture_one_circle pic(a, {0.2, 0.5}, false, square_cast_1, 100, 1e-10);
	draw_picture(pic, L"one_circle_4lab.bmp", pic_size);
	}
}