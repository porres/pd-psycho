
#include <math.h>
#include "m_pd.h"

static t_class *hz2bark_class;

typedef struct hz2bark{
    t_object  x_obj;
    t_outlet *x_outlet;
}t_hz2bark;

static void hz2bark_float(t_hz2bark *x, t_floatarg hz){
    float bark = 0;
    if(hz < 219.5)
        bark = 13.3 * atan(3*hz/4000); // Terhardt
    else{ // Traumuller
        bark = (26.81*hz)/(1960 + hz) - 0.53;
        if(bark > 20.1)
            bark = bark + 0.22 * (bark - 20.1);
    }
    outlet_float(x->x_outlet, bark);
}

static void *hz2bark_new(void){
    t_hz2bark *x = (t_hz2bark *)pd_new(hz2bark_class);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return (void *)x;
}

void hz2bark_setup(void){
    hz2bark_class = class_new(gensym("hz2bark"), (t_newmethod)hz2bark_new, 0,
                              sizeof(t_hz2bark), 0, 0);
	class_addfloat(hz2bark_class, hz2bark_float);
}
