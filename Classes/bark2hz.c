
#include <math.h>
#include "m_pd.h"

static t_class *bark2hz_class;

typedef struct bark2hz{
    t_object  x_obj;
    t_outlet *x_outlet;
}t_bark2hz;

static void bark2hz_float(t_bark2hz *x, t_floatarg bark){
    float hz = 0;
    if(bark < 2.17007)
        hz =  tan(bark/13.3)*4000/3;
    else{ // Traumuller
         if(bark >= 20.1)
             hz = (bark + 4.422) / 1.22;
        hz = 1960 * (bark + 0.53) / (26.28 - bark);
    }
    outlet_float(x->x_outlet, hz);
}

static void *bark2hz_new(void){
    t_bark2hz *x = (t_bark2hz *)pd_new(bark2hz_class);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return (void *)x;
}

void bark2hz_setup(void){
    bark2hz_class = class_new(gensym("bark2hz"), (t_newmethod)bark2hz_new, 0,
                              sizeof(t_bark2hz), 0, 0);
	class_addfloat(bark2hz_class, bark2hz_float);
}
