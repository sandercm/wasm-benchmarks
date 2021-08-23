#include <random>
#include <chrono>
#include <iostream>
#include <array>
#include <random>

#define MAXDIM 512
#define RUNS 10

using namespace std;

random_device rd;
mt19937 rng(rd());
uniform_int_distribution<int> uni(0, 100);

double randBetween() {
	return uni(rng);
}

class matrix
{
public:
	matrix(int dim)
		: dim_(dim)
	{}

	inline int dim() const {
		return dim_;
	}
	inline int& operator()(unsigned row, unsigned col) {
		return data_[dim_ * row + col];
	}

	inline int operator()(unsigned row, unsigned col) const {
		return data_[dim_ * row + col];
	}

private:
	int dim_;
	std::array<int, 9> data_{ {0, 0,0,0,0,0,0,0,0} };
};

void random_matrix_class(matrix& matrix) {
	for (int r = 0; r < matrix.dim(); r++)
		for (int c = 0; c < matrix.dim(); c++)
			matrix(r, c) = randBetween();
}

void mult_std(matrix const& a, matrix const& b, matrix& z) {
	for (int r = 0; r < a.dim(); r++) {
		for (int c = 0; c < a.dim(); c++) {
			z(r, c) = 0;
			for (int i = 0; i < a.dim(); i++)
				z(r, c) += a(r, i) * b(i, c);
		}
	}
}


int main(int n) {

	using chrono::high_resolution_clock;
	using chrono::duration_cast;
	using chrono::duration;
	using chrono::milliseconds;
	auto t1 = high_resolution_clock::now();

	matrix a(3);
	matrix b(3);
	matrix c(3);
	random_matrix_class(a);
	random_matrix_class(b);

	for (size_t i = 0; i < 100000; i++)
	{
		random_matrix_class(a);
		random_matrix_class(b);
		mult_std(a, b, c);
	}

	auto t2 = high_resolution_clock::now();

	duration<double, milli> ms_double = t2 - t1;
	cout << "duration was" << '\n';
	cout << ms_double.count() << '\n';
	return ms_double.count();
	return 0;
}