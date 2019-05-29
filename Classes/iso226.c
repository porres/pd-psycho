
#include <math.h>
#include "m_pd.h"

static t_class *iso226_class; /* not a t_func */

typedef struct iso226{
    t_object x_ob;
    t_outlet *x_outlet;
}t_iso226;

static void *iso226_new(void){
    t_iso226 *x = (t_iso226 *)pd_new(iso226_class);
    x->x_outlet = outlet_new(&x->x_ob, gensym("list"));
    return (void *)x;
}

static void iso226_float(t_iso226 *x, t_float ph){
    t_atom spl[29];
    int i;
    //f = [20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500];
    float af[] = {0.532, 0.506, 0.480, 0.455, 0.432, 0.409, 0.387, 0.367, 0.349, 0.330, 0.315,
                  0.301, 0.288, 0.276, 0.267, 0.259, 0.253, 0.250, 0.246, 0.244, 0.243, 0.243,
                  0.243, 0.242, 0.242, 0.245, 0.254, 0.271, 0.301};
    float Lu[] = {-31.6, -27.2, -23.0, -19.1, -15.9, -13.0, -10.3, -8.1, -6.2, -4.5, -3.1,
                   -2.0,  -1.1,  -0.4,   0.0,   0.3,   0.5,   0.0, -2.7, -4.1, -1.0,  1.7,
                    2.5,   1.2,  -2.1,  -7.1, -11.2, -10.7,  -3.1};
    float Tf[] = { 78.5,  68.7,  59.5,  51.1,  44.0,  37.5,  31.5,  26.5,  22.1,  17.9,  14.4,
                   11.4,   8.6,   6.2,   4.4,   3.0,   2.2,   2.4,   3.5,   1.7,  -1.3,  -4.2,
                   -6.0,  -5.4,  -1.5,   6.0,  12.6,  13.9,  12.3};
    if(ph < 0 || ph > 90){
        pd_error(x, "Phon value out of bounds!");
        return;
    }
    //Deriving sound pressure level from loudness level (iso226 sect 4.1)
    float Af[29];
    for(i = 0; i < 29; i++)
        Af[i] = 4.47E-3 * (pow(10,0.025*ph) - 1.15) + pow(0.4*pow(10,((Tf[i]+Lu[i])/10)-9 ), af[i]);
    for(i = 0; i < 29; i++)
        SETFLOAT(&spl[i], ((10/af[i]) * log10(Af[i])) - Lu[i] + 94);
    outlet_list(x->x_outlet,&s_list,29,spl);
}

static int eqlbandbins[43] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  11,  13,  15,  17,  19,
     22,  25,  28,  32,  36,  41,  46,  52,  58,  65,  73,  82,  92, 103,
    116, 129, 144, 161, 180, 201, 225, 251, 280, 312, 348, 388, 433, 483, 513
};

static float logfreqs[43];

void iso226_setup(void){
    iso226_class = class_new(gensym("iso226"), (t_newmethod)iso226_new, 0,
                             sizeof(t_iso226), 0, 0);
    class_addfloat(iso226_class,(t_method)iso226_float);
    for(int i = 0; i < 43; i++)
        logfreqs[i] = log(eqlbandbins[i]/512.0*44100.0); // ????????????
}
