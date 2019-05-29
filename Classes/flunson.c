
#include <math.h>
#include "m_pd.h"

//***** Convert Phons/hz to dB (Barlow's function, Robson/Dadson's data)
static float throfs(float hz){
    float khz, thr;
    khz = hz/1000;
    thr = 3.64*exp(-0.8*log(khz))+0.001 * pow(khz, 4);
    if(khz <= 15.41)
        thr = thr-6.5*exp(-0.6 * pow(khz-3.3, 2));
    return thr;
}

float e(float a){
    return fabs(a) < 80 ? exp(a) : 0;
}

float bell(float a, float b){
    return b != 0 ? e(-(a*a)/(b*16)) : 0;
}

float tanhyp(float x){
    float u = e(x);
    float v = e(-x);
    return (u+v) != 0 ? (u-v)/(u+v) : 0;
}

float hypt(float a, float b){
    return b != 0 ? 1-tanhyp(a/b) : 0;
}

float trait(int q, float phx){
    float rbf = 0;
    if (q == 1)  rbf = -711+573*hypt(phx+23,-137);
    if (q == 2)  rbf =  155+115*hypt(phx+23,62) + 648*hypt(phx-144,-35);
    if (q == 3)  rbf =   64+264*hypt(phx+23,114) + 53*hypt(phx-112,-32);
    if (q == 4)  rbf =   365+43*bell(phx-15,1419);
    if (q == 5)  rbf =  272+394*bell(phx-107,26) + 1241*hypt(phx+51,36);
    if (q == 6)  rbf =    26+28*bell(phx-108,44);
    if (q == 7)  rbf =   578+10*bell(phx-109,10);
    if (q == 8)  rbf =   86+103*bell(phx-109,22) + 61*hypt(phx+15,43);
    if (q == 9)  rbf =    62-34*bell(phx-63,66) + 39*hypt(phx-113,-162);
    if (q == 10) rbf =   -20+12*bell(phx-66,263) - 64*hypt(phx-153,-24);
    if (q == 11) rbf =    72-57*bell(phx-75,449) + 48*hypt(phx-116,-17);
    return rbf;
}

float dd(float h, float p){
    float x = 90.89 * log(fabs(h)) - 244;
    float phx = p * 1.03;
    return(0.288*(trait(1, phx) + trait(3, phx) * hypt(x-20, trait(2, phx))
                  + trait(6, phx) * bell(x-trait(4, phx), trait(5, phx))
                  + trait(9, phx) * bell(x-trait(7, phx), trait(8, phx))
                  + trait(11, phx) * hypt(x-620, trait(10, phx))));
}

// phon/hz to dbA
static float dbA(float hz, float ph){
    return dd (hz, ph) + (throfs(hz) - dd(hz,3)) * (hz/250000);
}

//**** Convert dBA and Hz to Phons (function by Clarence Barlow, data from Robson & Dadson)

static float phon(float hz, float db, int phonsteps){
/* NB: 1 dbA (the minimum workable value)
 gives -123.3 Ph @ 20 Hz, -30.0 Ph @ 20000 Hz */
	float ph, phbt, phtp, lphbt, lphtp; int i;
	float mb;
	phbt = -125.95;
    phtp = 150.00;
    i = 0;
	do{
		lphbt = phbt;
        lphtp = phtp;
		ph = (phbt+phtp)/2;
        mb = dbA(hz,ph);
		if(mb<db){
            phbt = ph;
            phtp = lphtp;
        }
		else{
            phbt = lphbt;
            phtp = ph;
        }
		i++;
	} while (i < phonsteps);
	return ph;
}

//****** [rmstodb] - Converts Linear Amplitude to dB(A)
static float gain2db(float amp){
    if(amp > 0.00001)
        return 100+20*log10f(amp);
    return 0;
}

//******* [dbtorms] - Converts dB(A) to Linear Amplitude
static float db2gain(float db){
    return pow(10, (db-100)/20);
}

//************************ FLUNSON OBJECT ************************/

t_class *flunson_class;

typedef struct flunson{
    t_object x_ob;
    t_outlet *x_outlet;
	int ac;
	float *amps;
}t_flunson;

void flunson_freqs(t_flunson *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(ac > x->ac)
        ac = x->ac; // too many args.
    t_atom outs[ac];
    for(int i = 0; i < ac; i++)
        SETFLOAT(outs+i, db2gain(phon(atom_getfloat(av+i), gain2db(x->amps[i]), 20)));
    outlet_list(x->x_outlet,&s_list,ac,outs);
}

void flunson_amps(t_flunson *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    // LATER see if resizebytes gets any faster
    if(x->amps)
        freebytes(x->amps, x->ac*sizeof(float));
    x->amps = (float *)getbytes(ac * sizeof(float));
	x->ac = ac;
    for(int i = 0; i < ac; i++)
		x->amps[i] = atom_getfloat(av+i);
}

void flunson_free (t_flunson *x){
    if(x->amps)
        freebytes(x->amps, x->ac*sizeof(float));
}

void *flunson_new(void){
    t_flunson *x = (t_flunson *)pd_new(flunson_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("list"), gensym("amps"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    return (void *)x;
}

void flunson_setup(void){
	flunson_class = class_new(gensym("flunson"), (t_newmethod)flunson_new,
                (t_method)flunson_free, sizeof(t_flunson), 0, 0);
	class_addlist(flunson_class, flunson_freqs);
	class_addanything(flunson_class, flunson_amps); // use class_addmethod with gensym("amps")
}
