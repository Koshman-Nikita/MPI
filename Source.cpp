#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[])
{
    double start, finish, duration;
    double h, x, lBound = 0.0, rBound = 1.0, IntervalLength = 0.0, currArea = 0.0, sum;

    int intervalNum = 0, procNum, procRank, namelen, i = 0, c = 0;

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);


    intervalNum = 200000000;

    if (procRank == 0)
    {
        cout << "Number of rectangles = " << intervalNum << endl;

        for (i = 0; i < procNum; ++i)
            printf("Process number: %d\n", i);
    }

    MPI_Bcast(&intervalNum, 1, MPI_INT, 0, MPI_COMM_WORLD);

    h = (rBound - lBound) / intervalNum;

    for (int i = procRank + 1; i <= intervalNum; i += procNum)
    {
        x = h * i;
        IntervalLength += sin(x + x * x * x);
        //c++;
    }
    //cout << c;

    currArea = h * IntervalLength;

    MPI_Reduce(&currArea, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    finish = MPI_Wtime();
    duration = finish - start;

    if (procRank == 0)
    {
        cout << std::setprecision(7) << "Result: " << sum << endl;
        cout << std::setprecision(7) << "Time of execution: " << duration << endl;
    }

    MPI_Finalize();

    return 0;
}