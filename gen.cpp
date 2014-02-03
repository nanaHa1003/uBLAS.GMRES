#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>
using namespace std;

int main(int argc, char **argv)
{
    int size = (argc == 2) ?atoi(argv[1]) :10 ;
    ofstream file("test.mtx");

    file << size << " " << size << " " << size * size  << std::endl;

    vector<double> base(size);
    #pragma omp parallel for
    for(int i = 0; i < base.size(); i++)
    {
        base[i] = (double) (i + 1);
    }

    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            file << i + 1 << " " << j + 1 << " ";
            file << pow(base[i], (double) j);
            file << std::endl;
        }
    }

    return 0;
}
