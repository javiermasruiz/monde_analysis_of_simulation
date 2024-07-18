
#include "TStopwatch.h"
#include "sim_analysis.h"

/*

Comentarios aqui

*/

void sim_analysis(){  

   TStopwatch t;

   readroot_to_txt();
   anger();
   plot_anger_root();

   t.Print();
}





