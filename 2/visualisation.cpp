#include "diagonal.h"
#include <windows.h>
#include <iomanip>

int main() {
	int n, m;
	int isOnlyLowTriangle;
	std::cout << "Matrix size: ";
	std::cin >> n;
	std::cout << "Diagonals count: ";
	std::cin >> m;
	std::cout << "Is only low triangle? (1/0) ";
	std::cin >> isOnlyLowTriangle;

	//-------------------------------------------------------------------------
	MatrixDiagonal a;
	generateDiagonalMatrix(n, 1, 10, generateRandomFormat(n, m), a);

	Matrix a_d;
	a.toDenseMatrix(a_d);

	matrix_diagonal_line_iterator mit(a.dimension(), a.getFormat(), isOnlyLowTriangle);

	auto drawMatrix = [&] () {
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
		system("cls");
		std::cout 
			<< "i = " << std::setw(3) << mit.i 
			<< "; j = " << std::setw(3) << mit.j 
			<< "; d = " << std::setw(3) << mit.d 
			<< "; dn = " << std::setw(3) << mit.dn 
			<< "; di = " << std::setw(3) << mit.di 
			<< std::endl
			<< std::endl;

		for (int i = 0; i < a_d.height(); ++i) {
			if (i == mit.i)
				std::cout << ">";
			for (int j = 0; j < a_d.width(); ++j) {
				if (i == mit.i && j == mit.j && !mit.isLineEnd())
					SetConsoleTextAttribute(hConsole, 12);
				else {
					SetConsoleTextAttribute(hConsole, 15);
					if (a_d(i, j) == 0)
						SetConsoleTextAttribute(hConsole, 8);
					if (i == j)
						SetConsoleTextAttribute(hConsole, 2);
				}

				if (i == mit.i && j == 0)
					std::cout << std::setw(1);
				else
					std::cout << std::setw(2);

				std::cout << a_d(i, j);
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;

		SetConsoleTextAttribute(hConsole, 15);
		if (!mit.isLineEnd())
			std::cout << " a[i, j] = " << a.begin(mit.dn)[mit.di] << std::endl;
		else
			std::cout << " Line end." << std::endl;

		system("pause");
	};

	//-------------------------------------------------------------------------
	for (; !mit.isEnd(); ++mit) {
		for (; !mit.isLineEnd(); ++mit)
			drawMatrix();
		drawMatrix();
	}
}