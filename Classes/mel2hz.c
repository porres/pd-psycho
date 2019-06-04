
#include <math.h>
#include "m_pd.h"

static t_class *mel2hz_class;

typedef struct mel2hz{
    t_object  x_obj;
    t_outlet *x_outlet;
}t_mel2hz;

static void mel2hz_float(t_mel2hz *x, t_floatarg mel){
    outlet_float(x->x_outlet, 700*(exp(mel/1127.01048)-1));
}

static void *mel2hz_new(void){
    t_mel2hz *x = (t_mel2hz *)pd_new(mel2hz_class);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return (void *)x;
}

void mel2hz_setup(void){
    mel2hz_class = class_new(gensym("mel2hz"), (t_newmethod)mel2hz_new,
            0, sizeof(t_mel2hz), 0, 0);
	class_addfloat(mel2hz_class, mel2hz_float);
}
