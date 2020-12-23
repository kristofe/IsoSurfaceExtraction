#include "exported_routines.h"
#include "IsoSurfaceExtraction.h"

void test(){
	printf("TESTING\n");
}

void run_quadratic_mc(const char* export_path, int res, float* voxelValues, float IsoValue, float* out_verts, int* out_vert_count, int* out_tris, int* out_tri_count)
{
    printf("running quadratic mc\n");
    extract_quadratic_isosurface(export_path, res, voxelValues, IsoValue, out_verts, out_vert_count, out_tris, out_tri_count);
}
