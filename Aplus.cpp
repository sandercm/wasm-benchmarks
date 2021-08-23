#include "math.h"
#include <array>
#include <chrono>
#include <cstring>
#include <queue>
#include <set>
#include <stack>
#include <tuple>
#include <iostream>
#include <emscripten/emscripten.h>
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
