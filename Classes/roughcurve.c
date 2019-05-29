
#include <math.h>
#include "m_pd.h"

#define EXP_1 exp(1)

static t_class *roughcurve_class;

typedef struct roughcurve{
    t_object  x_obj;
    t_outlet *x_outlet;
    int       x_curve;
}t_roughcurve;

static float parncutt(float bark_dif){ // Curve = 0
    //    return bark_dif > 1.2 ? 0 : pow(4*EXP_1*bark_dif / exp(4*bark_dif), 2);
    return pow(4*EXP_1*bark_dif / exp(4*bark_dif), 2);
}

static float sethares(float bark_dif){ // Curve = 1
    return 5.56309 * (exp(bark_dif * -3.51) - exp(bark_dif * -5.75));
}

static void roughcurve_float(t_roughcurve *x, t_floatarg f){
    outlet_float(x->x_outlet, x->x_curve ? sethares(f) : parncutt(f));
}

static void roughcurve_curve(t_roughcurve *x, t_float f){
    x->x_curve = (int)(f != 0);
}

static void *roughcurve_new(t_floatarg f){
    t_roughcurve *x = (t_roughcurve *)pd_new(roughcurve_class);
    x->x_curve = (int)(f != 0);
    x->x_outlet = outlet_new(&x->x_obj, 0);
    return(void *)x;
}

void roughcurve_setup(void){
    roughcurve_class = class_new(gensym("roughcurve"), (t_newmethod)roughcurve_new, 0,
                              sizeof(t_roughcurve), 0, A_DEFFLOAT, 0);
	class_addfloat(roughcurve_class, roughcurve_float);
    class_addmethod(roughcurve_class, (t_method)roughcurve_curve,
                    gensym("curve"), A_FLOAT, 0);
}
