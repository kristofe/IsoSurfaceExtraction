#if 1
//This is just so cffi doesn't error on parsing this file
#include <stdint.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

void run_quadratic_mc(const char* export_path, int res, float* voxelValues, float IsoValue, float* out_verts, int* out_vert_count, int* out_tris, int* out_tri_count);
void test();

#ifdef __cplusplus
}
#endif
