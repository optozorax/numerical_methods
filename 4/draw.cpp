#include <sstream>
#include <fstream>
#include <iomanip>
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
	void init_solve(const pair<sle_f, fnmv_f>& p, const xn_t& x_0, bool is_numeric, const sqr_f& square_cast, int iterations, double residual) {
		m_x_0 = x_0;
		m_iterations = iterations;
		m_residual = residual;
		if (is_numeric) {
			auto j = calc_jacobi_matrix_numeric_functon(p.second);
			auto sle = get_sle_function(j, p.second);
			sle = square_cast(sle);
			m_solved = solve(sle, x_0, iterations, residual, false);
			m_f = get_f(sle);
		} else {
			auto sle = square_cast(p.first);
			m_solved = solve(sle, x_0, iterations, residual, false);
			m_f = get_f(sle);
		}
	}
	virtual void init_brd(FindBorders& brd) const {
		brd.process(vec2(0, 0));
		for (int i = 0; i < m_solved.x_process.size(); ++i)
			brd.process(vec2(m_solved.x_process[i][0], m_solved.x_process[i][1]));
	}
	virtual double getResidual(vec2 a) const {
		return length(m_f({ a.x, a.y }));
	}
	virtual void draw(ImageDrawing_aa& img, const FindBorders& brd) const = 0;
	vector<vec2> getProcess(void) const {
		vector<vec2> result;
		for (auto& i : m_solved.x_process)
			result.push_back(vec2(i[0], i[1]));
		return result;
	}

	xn_t get_x_0(void) const { return m_x_0; }
	int get_iterations(void) const { return m_iterations; }
	double get_residual(void) const { return m_residual; }
	solved_t get_solved(void) const { return m_solved; }
protected:
	solved_t m_solved;
	fnm_f m_f;

	xn_t m_x_0;
	int m_iterations;
	double m_residual;
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
		return length(m_f({a.x, a.y, r}));
	}

	void draw(ImageDrawing_aa& img, const FindBorders& brd) const {
		img.drawPolyline(calc_circle(m_a, brd));
		img.drawPolyline(calc_circle(m_b, brd));
		img.drawPolyline(calc_circle(m_c, brd));

		Pen oldPen = img.getPen();
		oldPen.clr = Orange;
		img.setPen(oldPen);
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
			pos = sqrt(pos);
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

	img.setPen(Pen(2, Blue));
	pic.draw(img, brd);

	// Рисуем линии поиска решения
	auto process = pic.getProcess();
	for (int i = 0; i < process.size()-1; i++) {
		auto a = brd.toImg(process[i]);
		auto b = brd.toImg(process[i+1]);

		img.setPen(Pen(3.5, setAlpha(Green, 100)));
		img.drawLine(a, b);
	}

	// Рисуем точки поиска решения
	for (int i = 0; i < process.size(); i++) {
		auto a = brd.toImg(process[i]);

		img.setBrush(Brush(setAlpha(twg::Miku, 200)));
		img.drawPolygon(computeEllipse(Point_d(4, 4)).move(a));
	}

	saveToBmp(&img, filename + L".bmp");

	// TODO добавить сюда вывод в файл процесса
	ofstream fout(filename + L".dat");

	auto solved = pic.get_solved();
	fout << "k\tbeta\tresidual\t";
	if (solved.point.size() == 2)
		fout << "x\ty" << endl;
	else
		fout << "x\ty\tr" << endl;
	for (int i = 1; i < solved.x_process.size(); i++) {
		fout << i << setprecision(2) << scientific << "\t" << solved.beta_process[i] << "\t" << solved.residual_process[i] << "\t" << setprecision(10) << fixed; 
		auto point = solved.x_process[i];
		for (int j = 0; j < point.size() - 1; j++)
			fout << point[j] << "\t";
		fout << point.back() << endl;
	}

	fout << setprecision(10) << fixed;
	fout << "# x_0 = " << pic.get_x_0() << endl;
	fout << "# max_iter = " << pic.get_iterations() << endl;
	fout << setprecision(2) << scientific;
	fout << "# required residual = " << pic.get_residual() << endl;

	fout << endl << "# ";

	switch (solved.exit_type) {
		case EXIT_ITER: fout << "Exit by iterations" << endl; break;
		case EXIT_RESIDUAL: fout << "Exit by residual" << endl; break;
		case EXIT_STEP: fout << "Exit by step" << endl; break;
		case EXIT_BETA: fout << "Exit by beta" << endl; break; 
		case EXIT_ERROR: fout << "Exit by error" << endl; break;
	}

	fout << "# iter = " << pic.get_solved().iterations << endl;
	fout << setprecision(2) << scientific;
	fout << "# Residual = " << pic.get_solved().residual << endl;
	fout << setprecision(10) << fixed;
	fout << "# Result point = " << pic.get_solved().point << endl;

	fout.close();
}

int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow) {
	const double required_residual = 1e-10;
	const int maxiter = 30;
	#ifdef _DEBUG
		int pic_size = 50;
	#else
		int pic_size = 700;
	#endif

	// Делаем таблицы для сходимости метода в зависимости от размера шага при взятии производной
	for (int i = 0; i < 10; i++) {
		partial_derivative_step = pow(2.0, double(-i));
		circle a(0, 0, 0.8), b(sqrt(2.0) / 2, sqrt(2.0) / 2, 0.2);
		line c = { point(sqrt(2.0)*(1.0 / 2.0 - 0.1), sqrt(2.0)*(1.0 / 2.0 - 0.1)), point(-1, 0) };
		Picture_two_circles_and_line pic(a, b, c, { -1, 0.5 }, true, square_cast_3, maxiter, required_residual);
		draw_picture(pic, L"img/two_circles_and_line_3_step_" + to_wstring(i) + L"_4lab_numeric", pic_size);
	}

	partial_derivative_step = 1e-9;

	auto make_images = [&](bool is_numeric, wstring append) {
		/*{
			circle a(0, 0, 1), b(5, 0, 1), c(10, 5, 5);
			Picture_three_circles pic(a, b, c, false, false, false, { -1, -1, 13 }, is_numeric, square_cast_none, 100, 1e-10);
			draw_picture(pic, L"img/three_circles_4lab" + append, pic_size);
		}*/

		auto draw_three_lines = [&](int no, const sqr_f& sq) {
			vec2 a(0, 0), b(1, 2), c(3, 1);
			Picture_three_lines pic(a, b, c, { -3, -3 }, is_numeric, sq, maxiter, required_residual);
			draw_picture(pic, L"img/three_lines_4lab_" + to_wstring(no) + append, pic_size);
		};

		draw_three_lines(2, square_cast_2);
		draw_three_lines(3, square_cast_3);
		draw_three_lines(4, square_cast_4);

		{
			vec2 a(0.5, 0.5), b(3, 1), c(-1, 3);
			Picture_three_lines pic(a, b, c, { -3, -3 }, is_numeric, composition(square_cast_3, square_cast_mul({1, 1, 10})), 100, required_residual);
			draw_picture(pic, L"img/three_lines_4lab_3_weight" + append, pic_size);
		}

		{
			vec2 a(2, 2);
			Picture_sin_and_line pic(a, { -1, 1 }, is_numeric, square_cast_none, maxiter, required_residual);
			draw_picture(pic, L"img/sin_and_line_4lab" + append, pic_size);
		}

		{
			circle a(0, 0, 1), b(1, 0.7, 1.2); // Пересекаются в двух точках
			Picture_two_circles pic(a, b, { -1, 0.5 }, is_numeric, square_cast_none, maxiter, required_residual);
			draw_picture(pic, L"img/two_circles_1_4lab" + append, pic_size);
		}

		{
			circle a(0, 0, 0.8), b(sqrt(2.0) / 2, sqrt(2.0) / 2, 0.2); // Пересекаются в одной точке
			{
				// Начальное приближение лежит на оси симметрии
				Picture_two_circles pic(a, b, { 0.156405, 0.974966 }, is_numeric, square_cast_none, maxiter, required_residual);
				draw_picture(pic, L"img/two_circles_2_1_4lab" + append, pic_size);
			}
			{
				// Начальное приближение лежит на линии, соединяющей центры окружностей
				Picture_two_circles pic(a, b, { -sqrt(2.0), -sqrt(2.0)+0.01 }, is_numeric, square_cast_none, maxiter, required_residual);
				draw_picture(pic, L"img/two_circles_2_2_4lab" + append, pic_size);
			}
			{
				// Начальное приближение лежит в центре одной из окружностей
				Picture_two_circles pic(a, b, { b.c.x, b.c.y + 0.01 }, is_numeric, square_cast_none, maxiter, required_residual);
				draw_picture(pic, L"img/two_circles_2_3_4lab" + append, pic_size);
			}
		}

		{
			circle a(0, 0, 0.8), b(1, 0.7, 0.3); // Не пересекаются
			Picture_two_circles pic(a, b, { -1, 0.5 }, is_numeric, square_cast_none, maxiter, required_residual);
			draw_picture(pic, L"img/two_circles_3_4lab" + append, pic_size);
		}

		auto draw_two_circles_and_line = [&](int no, const sqr_f& sq) {
			circle a(0, 0, 0.8), b(sqrt(2.0) / 2, sqrt(2.0) / 2, 0.2);
			line c = { point(sqrt(2.0)*(1.0 / 2.0 - 0.1), sqrt(2.0)*(1.0 / 2.0 - 0.1)), point(-1, 0) };
			Picture_two_circles_and_line pic(a, b, c, { -1, 0.5 }, is_numeric, sq, maxiter, required_residual);
			draw_picture(pic, L"img/two_circles_and_line_" + to_wstring(no) + L"_4lab" + append, pic_size);
		};

		draw_two_circles_and_line(2, square_cast_2);
		draw_two_circles_and_line(3, square_cast_3);
		draw_two_circles_and_line(4, square_cast_4);

		{
			circle a(0, 0, 1);
			Picture_one_circle pic(a, { 0.2, 0.5 }, is_numeric, square_cast_1, maxiter, required_residual);
			draw_picture(pic, L"img/one_circle_4lab_1" + append, pic_size);
		}
	};
	make_images(false, L"_analytic");
	make_images(true, L"_numeric");
}