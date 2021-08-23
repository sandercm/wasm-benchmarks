#include <math.h>
#include <iostream>
#include <chrono>
#include <array>
#include <random>
#include <queue>
#include <set>
#include <stack>
#include <tuple>
#include <iostream>
#include <emscripten/emscripten.h>
extern "C"
{

    struct Sphere
    {
        double x, y, z, radius;
        Sphere(double _x, double _y, double _z, double _radius)
        {
            x = _x;
            y = _y;
            z = _z;
            radius = _radius;
        }
    };
    std::random_device rd;                               // only used once to initialise (seed) engine
    std::mt19937 rng(rd());                              // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uni(-1000, 1000); // guaranteed unbiased
    double randBetween()
    {
        return uni(rng);
    }

    bool sphereToSphereCollision(struct Sphere sphere1, struct Sphere sphere2)
    {
        return (sphere1.radius + sphere2.radius) > hypot(hypot(sphere1.x - sphere2.x, sphere1.y - sphere2.y), sphere1.z - sphere2.z);
    }

    EMSCRIPTEN_KEEPALIVE
    double runSphereToSphereCollision(int loops)
    {
        using std::chrono::duration;
        using std::chrono::duration_cast;
        using std::chrono::high_resolution_clock;
        using std::chrono::milliseconds;

        std::vector<Sphere> spheres;
        for (size_t i = 0; i < loops * 2; i++)
        {
            spheres.push_back({randBetween(), randBetween(), randBetween(), randBetween()});
        }
        auto t1 = high_resolution_clock::now();
        double total = 0.0;
        for (size_t i = 0; i < loops; i += 2)
        {
            total += sphereToSphereCollision(spheres[i], spheres[i + 1]);
        }
        auto t2 = high_resolution_clock::now();

        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = t2 - t1;
        std::cout << "duration was ";
        std::cout << ms_double.count() << '\n';
        std::cout << "total was ";
        std::cout << total << '\n';
        return ms_double.count();
    }
    // ----------------------
    // QUICKSORT
    // ----------------------

    void swap(int *a, int *b)
    {
        int t = *a;
        *a = *b;
        *b = t;
    }

    int partition(int arr[], int low, int high)
    {
        int pivot = arr[high]; // pivot
        int i = (low - 1);     // Index of smaller element

        for (int j = low; j <= high - 1; j++)
        {
            // If current element is smaller than or
            // equal to pivot
            if (arr[j] <= pivot)
            {
                i++; // increment index of smaller element
                swap(&arr[i], &arr[j]);
            }
        }
        swap(&arr[i + 1], &arr[high]);
        return (i + 1);
    }

    void quickSort(int arr[], int low, int high)
    {
        if (low < high)
        {
            int pi = partition(arr, low, high);

            quickSort(arr, low, pi - 1);
            quickSort(arr, pi + 1, high);
        }
    }

    EMSCRIPTEN_KEEPALIVE
    double runQuickSort(int n, int size)
    {
        using std::chrono::duration;
        using std::chrono::duration_cast;
        using std::chrono::high_resolution_clock;
        using std::chrono::milliseconds;

        auto t1 = high_resolution_clock::now();
        std::vector<int> vec(size);
        for (size_t i = 0; i < n; i++)
        {
            for (int x = 0; x < size; ++x)
                vec[x] = randBetween();
            quickSort(&vec[0], 0, vec.size() - 1);
        }

        auto t2 = high_resolution_clock::now();

        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = t2 - t1;
        std::cout << "duration was" << '\n';
        std::cout << ms_double.count() << '\n';
        return ms_double.count();
    }
}

// ----------------------
// A* algorithm
// ----------------------

using namespace std;


typedef pair<int, int> Pair;
typedef tuple<double, int, int> Tuple;


struct cell {

	Pair parent;
	double f, g, h;

	cell()
		: parent(-1, -1)
		, f(-1)
		, g(-1)
		, h(-1)
	{
	}
};


template <size_t ROW, size_t COL>
bool isValid(const array<array<int, COL>, ROW>& grid,
	const Pair& point)
{
	if (ROW > 0 && COL > 0)
		return (point.first >= 0) && (point.first < ROW)
		&& (point.second >= 0)
		&& (point.second < COL);

	return false;
}


template <size_t ROW, size_t COL>
bool isUnBlocked(const array<array<int, COL>, ROW>& grid,
	const Pair& point)
{
	return isValid(grid, point)
		&& grid[point.first][point.second] == 1;
}

bool isDestination(const Pair& position, const Pair& dest)
{
	return position == dest;
}

double calculateHValue(const Pair& src, const Pair& dest)
{
	return sqrt(pow((src.first - dest.first), 2.0)
		+ pow((src.second - dest.second), 2.0));
}

template <size_t ROW, size_t COL>
void aStarSearch(const array<array<int, COL>, ROW>& grid,
	const Pair& src, const Pair& dest)
{
	if (!isValid(grid, src)) {
		printf("Source is invalid\n");
		return;
	}

	if (!isValid(grid, dest)) {
		printf("Destination is invalid\n");
		return;
	}

	if (!isUnBlocked(grid, src)
		|| !isUnBlocked(grid, dest)) {
		printf("Source or the destination is blocked\n");
		return;
	}

	if (isDestination(src, dest)) {
		printf("We are already at the destination\n");
		return;
	}

	bool closedList[ROW][COL];
	memset(closedList, false, sizeof(closedList));

	array<array<cell, COL>, ROW> cellDetails;

	int i = src.first;
	int	j = src.second;

	cellDetails[i][j].f = 0.0;
	cellDetails[i][j].g = 0.0;
	cellDetails[i][j].h = 0.0;
	cellDetails[i][j].parent = { i, j };

	std::priority_queue<Tuple, std::vector<Tuple>, std::greater<Tuple>> openList;

	openList.emplace(0.0, i, j);

	while (!openList.empty()) {
		std::tie(std::ignore, i, j) = openList.top();

		openList.pop();
		closedList[i][j] = true;

		for (int add_x = -1; add_x <= 1; add_x++) {
			for (int add_y = -1; add_y <= 1; add_y++) {
				Pair neighbour(i + add_x, j + add_y);

				if (!isValid(grid, neighbour)) {
					return;
				}

				if (isDestination(neighbour, dest)) {
					cellDetails[neighbour.first][neighbour.second].parent = { i, j };
					return;
				}

				else if (!closedList[neighbour.first][neighbour.second] && isUnBlocked(grid, neighbour)) {
					double gNew = cellDetails[i][j].g + 1.0;
					double hNew = calculateHValue(neighbour, dest);
					double fNew = gNew + hNew;

					if (cellDetails[neighbour.first][neighbour.second].f == -1 || cellDetails[neighbour.first][neighbour.second].f > fNew) {
						openList.emplace(fNew, neighbour.first, neighbour.second);
						cellDetails[neighbour.first][neighbour.second].g = gNew;
						cellDetails[neighbour.first][neighbour.second].h = hNew;
						cellDetails[neighbour.first][neighbour.second].f = fNew;
						cellDetails[neighbour.first][neighbour.second].parent = { i, j };
					}
				}
			}
		}
	}
	printf("Failed to find the Destination Cell\n");
}
extern "C" {
	// Driver program to test above function
	EMSCRIPTEN_KEEPALIVE
	double runAStar(int n)
	{
		array<array<int, 10>, 9> grid{
			{ { { 1, 0, 1, 1, 1, 1, 0, 1, 1, 1 } },
			  { { 1, 1, 1, 0, 1, 1, 1, 0, 1, 1 } },
			  { { 1, 1, 1, 0, 1, 1, 0, 1, 0, 1 } },
			  { { 0, 0, 1, 0, 1, 0, 0, 0, 0, 1 } },
			  { { 1, 1, 1, 0, 1, 1, 1, 0, 1, 0 } },
			  { { 1, 0, 1, 1, 1, 1, 0, 1, 0, 0 } },
			  { { 1, 0, 0, 0, 0, 1, 0, 0, 0, 1 } },
			  { { 1, 0, 1, 1, 1, 1, 0, 1, 1, 1 } },
			  { { 1, 1, 1, 0, 0, 0, 1, 0, 0, 1 } } }
		};

		Pair src(8, 0);
		Pair dest(0, 0);

		using std::chrono::high_resolution_clock;
		using std::chrono::duration_cast;
		using std::chrono::duration;
		using std::chrono::milliseconds;
		auto t1 = high_resolution_clock::now();
		for (size_t i = 0; i < n; i++)
		{
			aStarSearch(grid, src, dest);
		}
		auto t2 = high_resolution_clock::now();

		/* Getting number of milliseconds as a double. */
		duration<double, std::milli> ms_double = t2 - t1;
		std::cout << "duration was ";
		std::cout << ms_double.count() << '\n';
		return ms_double.count();
	}
}

// --------------------------------
// Matrix algo
// --------------------------------
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

extern "C" {
	EMSCRIPTEN_KEEPALIVE
	double runMatrixMult(int n) {

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

	for (size_t i = 0; i < n; i++)
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
}
}