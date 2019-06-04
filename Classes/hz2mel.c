
#include <math.h>
#include "m_pd.h"

static t_class *hz2mel_class;

typedef struct hz2mel{
    t_object  x_obj;
    t_outlet *x_outlet;
}t_hz2mel;

static void hz2mel_float(t_hz2mel *x, t_floatarg hz){
    outlet_float(x->x_outlet, 1127.01048*log(1+(hz/700)));
}

static void *hz2mel_new(void){
    t_hz2mel *x = (t_hz2mel *)pd_new(hz2mel_class);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return (void *)x;
}

void hz2mel_setup(void){
    hz2mel_class = class_new(gensym("hz2mel"), (t_newmethod)hz2mel_new,
            0, sizeof(t_hz2mel), 0, 0);
	class_addfloat(hz2mel_class, hz2mel_float);
}
