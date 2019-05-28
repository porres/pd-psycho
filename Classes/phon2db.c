
#include <string.h>
#include <math.h>
#include "m_pd.h"

// functions
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

static float db(float hz, float ph){ // convert to db
    return dd(hz, ph) + (throfs(hz) - dd(hz,3)) * (hz/250000);
}

static int eqlbandbins[43]= {
    1,  2,  3,  4,  5,  6,  7,  8,  9, 11, 13, 15, 17, 19,
    22, 25, 28, 32, 36, 41, 46, 52, 58, 65, 73, 82, 92,103,
    116,129,144,161,180,201,225,251,280,312,348,388,433,483,513
};

static float logfreqs[43];

// object
typedef struct phon2db{
    t_object    x_obj;
    t_outlet   *x_outlet;
    float       x_hz;
}t_phon2db;

static t_class *phon2db_class;

static void *phon2db_new(t_float hz){
    t_phon2db *x = (t_phon2db *)pd_new(phon2db_class);
    x->x_hz = hz;
    floatinlet_new(&x->x_obj, &x->x_hz);
    x->x_outlet = outlet_new(&x->x_obj, gensym("float"));
    return (void *)x;
}

static void phon2db_float(t_phon2db *x, float phon){
    outlet_float(x->x_outlet, db(x->x_hz, phon));
}

void phon2db_setup(void){
    phon2db_class = class_new(gensym("phon2db"), (t_newmethod)phon2db_new, 0,
                    sizeof(t_phon2db), 0, A_DEFFLOAT, 0);
    class_addfloat(phon2db_class, (t_method)phon2db_float);
    for(int i = 0; i < 43; i++)
        logfreqs[i] = log(eqlbandbins[i]/512.0*44100.0); // ????????????
}
