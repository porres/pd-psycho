
#include <m_pd.h>
#include <math.h>

typedef struct{
	t_object x_obj;
}t_indigestibility;

static t_class *indigestibility_class;

void factorise(int n, int *acp, t_atom *av){
    *acp = 0;
    int p = 2;
    while(n > 1){
        if(n % p == 0){
            SETFLOAT(&av[*acp], p);
            (*acp)++;
            n /= p;
        }
        else
        p++;
        if(*acp >= 24)
        break; // safeguard
    }
}

float indigestibility(int n){
    int ac = 0;
	t_atom av[24];
	factorise(n, &ac, av);
	float ind = 0;
	for(int i = 0; i < ac; i++){
		float p = atom_getfloat(&av[i]);
		ind += 2 * (p-1) * (p-1)/p;
	}
	return ind;
}

void indigestibility_float (t_indigestibility *x, t_floatarg f){
	if (f <= 0 || f > 16777216){
		pd_error(x, "number %f out of range (0 ... 16777216)", f);
		return;
	}
    outlet_float(x->x_obj.te_outlet, indigestibility(f));
}

void *indigestibility_new(void){
	t_indigestibility *x = (t_indigestibility *)pd_new(indigestibility_class);
	outlet_new((t_object *)x, 0);
	return(x);
}

void indigestibility_setup (void){    
	indigestibility_class = class_new(gensym("indigestibility"), (t_newmethod)indigestibility_new, 0,
        sizeof(t_indigestibility), 0, 0);
	class_addfloat(indigestibility_class,indigestibility_float);
}
