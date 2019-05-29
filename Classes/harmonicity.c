
#include <m_pd.h>
#include <math.h>

typedef struct{
	t_object x_obj;
	t_float x_g;
}t_harmonicity;

static t_class *harmonicity_class;

void factorise(int n, int *argcp, t_atom *av){
    *argcp = 0;
    int p = 2;
    while(n > 1){
        if(n % p == 0){
            SETFLOAT(&av[*argcp], p);
            (*argcp)++;
            n /= p;
        }
        else
            p++;
        if(*argcp >= 24)
            break; // safeguard
    }
}

float indigestibility(int n){
    int argc = 0;
    t_atom av[24];
    factorise(n, &argc, av);
    int i;
    float ind = 0;
    for(i = 0; i < argc; i++){
        float p = atom_getfloat(&av[i]);
        ind += 2 * (p-1)*(p-1)/p;
    }
    return ind;
}

int sgn(int v){
    return v < 0 ? -1 : v > 0 ? 1 : 0;
}

void harmonicity_float (t_harmonicity *x, t_floatarg f){
	float g = x->x_g;
	if(f <= 0 || f > 16777216){
        pd_error(x, "left inlet number %f out of range (1 ... 16777216)", f);
        return;
    }
	if(g <= 0 || g > 16777216){
        pd_error(x, "left inlet number %f out of range (1 ... 16777216)", g);
        return;
    }
    outlet_float(x->x_obj.te_outlet, sgn(g-f)/(indigestibility(g)+indigestibility(f)));
}

void *harmonicity_new(t_float g){
	t_harmonicity *x = (t_harmonicity *)pd_new(harmonicity_class);
    x->x_g = g;
    floatinlet_new(&x->x_obj, &x->x_g);
	outlet_new((t_object *)x, 0);
	return(x);
}

void harmonicity_setup (void){
	harmonicity_class = class_new(gensym("harmonicity"), (t_newmethod)harmonicity_new, 0,
            sizeof(t_harmonicity), 0, 0);
	class_addfloat(harmonicity_class,harmonicity_float);
}
