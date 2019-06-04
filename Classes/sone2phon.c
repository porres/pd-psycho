
#include <math.h>
#include "m_pd.h"

#define lOG10_2 log10(2)

static t_class *sone2phon_class;

typedef struct sone2phon{
    t_object  x_obj;
    t_outlet *x_outlet;
}t_sone2phon;

static void sone2phon_float(t_sone2phon *x, t_floatarg ph){
    outlet_float(x->x_outlet, 40 + 33.22*log10f(ph));
}

static void *sone2phon_new(void){
    t_sone2phon *x = (t_sone2phon *)pd_new(sone2phon_class);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return (void *)x;
}

void sone2phon_setup(void){
    sone2phon_class = class_new(gensym("sone2phon"), (t_newmethod)sone2phon_new, 0,
                              sizeof(t_sone2phon), 0, 0);
	class_addfloat(sone2phon_class, sone2phon_float);
}
