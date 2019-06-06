// by porres, 2019, based on "Harmony: A Psychoacoustical Approach" (Parncutt, R. 1989).

#include "m_pd.h"
#include <math.h>

typedef struct floatlist{
	int         n;
	t_float    *v;
}t_floatlist;

typedef struct yl{
    t_object    x_obj;
    t_floatlist amps;
    t_outlet   *outlet0;
    t_outlet   *outlet1;
}t_yl;

t_class *yl_class;

void yl_list(t_yl *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(ac != x->amps.n){
        pd_error(x, "[yl]: input lists don't have the same length");
        return;
    }
    t_atom at[ac];
    for(int i = 0; i < ac; i++){
        float kHz = atom_getfloat(av+i) * 0.001;
        float yl = 3.64 * pow(kHz, -0.8) - 6.5 * exp(-0.6 * pow(kHz -3.3, 2))
        + 0.001 * pow(kHz, 4); // Eq.2
        SETFLOAT(&at[i], rmstodb(x->amps.v[i])-yl); // YL = dB-yl
    }
    outlet_list(x->outlet1, &s_list, ac, at);
    outlet_list(x->outlet0, &s_list, ac, av);
}

void yl_amps(t_yl *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->amps.v)
        freebytes(x->amps.v, x->amps.n*sizeof(float));
    x->amps.v = (float *)getbytes(ac*sizeof(float));
    x->amps.n = ac;
    for(int i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(av+i);
}

void yl_free(t_yl *x){
    freebytes(x->amps.v, x->amps.n*sizeof(float));
}

void *yl_new(void){
    t_yl *x = (t_yl *)pd_new(yl_class);
    x->amps.n = 0;
    x->amps.v = 0;
    inlet_new((t_object *)x, (t_pd *)x, gensym("list"), gensym("amps"));
    x->outlet0 = outlet_new((t_object *)x, gensym("list"));
    x->outlet1 = outlet_new((t_object *)x, gensym("list"));
	return(x);
}

void yl_setup(void){ // (Hz, linamp) -> (Hz, YL)
    yl_class = class_new(gensym("yl"), (t_newmethod)yl_new, (t_method)yl_free,
            sizeof(t_yl), 0, 0);
	class_addlist(yl_class, (t_method)yl_list);
	class_addmethod(yl_class, (t_method)yl_amps, gensym("amps"), A_GIMME, 0);
}
