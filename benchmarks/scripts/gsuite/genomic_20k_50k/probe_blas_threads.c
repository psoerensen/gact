/* Read thread controls from the actual loaded Windows DLLs; no package changes. */
#include <windows.h>
__declspec(dllexport) void probe_blas_threads(int *counts, int *parallel) {
 const char *names[]={"Rblas.dll","Rlapack.dll","libopenblas.dll"};
 for(int i=0;i<3;i++){
  HMODULE h=GetModuleHandleA(names[i]);
  int (*get_threads)(void)=h?(int (*)(void))GetProcAddress(h,"openblas_get_num_threads"):0;
  int (*get_parallel)(void)=h?(int (*)(void))GetProcAddress(h,"openblas_get_parallel"):0;
  counts[i]=get_threads?get_threads():-1;
  parallel[i]=get_parallel?get_parallel():-1;
 }
}
