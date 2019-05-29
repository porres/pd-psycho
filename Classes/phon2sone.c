
#include <math.h>
#include "m_pd.h"

#define lOG10_2 log10(2)

static t_class *phon2sone_class;

typedef struct phon2sone{
    t_object  x_obj;
    t_outlet *x_outlet;
}t_phon2sone;

static void phon2sone_float(t_phon2sone *x, t_floatarg ph){
    float sone = 0;
    if(ph >  40)
        sone = pow(10, lOG10_2*(ph-40)/10);
    else if(ph > 1)
        sone = pow(ph/40, 2.86);
    outlet_float(x->x_outlet, sone);
}

static void *phon2sone_new(void){
    t_phon2sone *x = (t_phon2sone *)pd_new(phon2sone_class);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return (void *)x;
}

void phon2sone_setup(void){
    phon2sone_class = class_new(gensym("phon2sone"), (t_newmethod)phon2sone_new, 0,
                              sizeof(t_phon2sone), 0, 0);
	class_addfloat(phon2sone_class, phon2sone_float);
}
