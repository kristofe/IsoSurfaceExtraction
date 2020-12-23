#include "exported_routines.h"
#include "IsoSurfaceExtraction.h"

void test(){
	printf("TESTING\n");
}

void run_quadratic_mc(const char* export_path, int res, float* voxelValues, float IsoValue)
{
    printf("running quadratic mc\n");
    extract_quadratic_isosurface(export_path, res, voxelValues, IsoValue);
}
