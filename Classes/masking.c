// by porres, 2019, based on "Harmony: A Psychoacoustical Approach" (Parncutt, R. 1989).

#include "m_pd.h"
#include <math.h>

typedef struct floatlist{
	int n;
	t_float *v;
}t_floatlist;

typedef struct masking{
    t_object    x_obj;
    t_floatlist amps;
    t_outlet   *outlet0;
    t_outlet   *outlet1;
    float       x_kM;
}t_masking;

t_class *masking_class;

void masking_kM(t_masking *x, float f){
    x->x_kM = f;
}

float erb(float Hz){
    float kHz = Hz/1000;
    return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

void masking_list(t_masking *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(ac != x->amps.n){
        pd_error(x, "[masking]: input lists don't have the same length");
        return;
    }
    int i, j;
    for(i = 0; i < ac; i++){
        float kHz = atom_getfloat(av+i) * 0.001;
        float LTh = 3.64 * pow(kHz, -0.8) - 6.5 * exp(-0.6 * pow(kHz -3.3, 2))
        + 0.001 * pow(kHz, 4);  // Eq.2
        float dB = rmstodb(x->amps.v[i]);
        x->amps.v[i] = dB-LTh; // YL
    }
    t_atom at[ac];
    for(i = 0; i < ac; i++){
        float YL_i = x->amps.v[i];
        float sum = 0;
        for(j = 0; j < ac; j++){
            if(i != j){ // masker (which doesn't mask itself)
                float YL_j = x->amps.v[j];
                float dif = fabs(erb(atom_getfloat(av+j)) - erb(atom_getfloat(av+i)));
                // pitch difference in critical bandwidths
                if(dif < 3.){ // otherwise masking is negligible
                    // (Eq.4) Partial Masking Level (masking due to one masker)
                    float PML = YL_j - x->x_kM * dif;
                    // The masking gradient kM is typically about 12 dB per cb
                    // Assumption: "Typical" complex tones in reg 4 have 10 audible harmonics.
                    sum += pow(10, PML/20.); // add amplitudes
                }
            }
        }
        // (Eq.5) Overall Masking Level (due to all maskers)
        float ML = fmax(20.*log10(sum), 0);
        // (Eq.6) Audible Level (level above masked threshold)
        float AL = fmax(YL_i-ML, 0);
        SETFLOAT(&at[i], AL);
    }
    outlet_list(x->outlet1, &s_list, ac, at);
    outlet_list(x->outlet0, &s_list, ac, av);
}

void masking_amps(t_masking *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->amps.v)
        freebytes(x->amps.v, x->amps.n*sizeof(float));
    x->amps.v = (float *)getbytes(ac*sizeof(float));
    x->amps.n = ac;
    for(int i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(av+i);
}

void *masking_new(void){
    t_masking *x = (t_masking *)pd_new(masking_class);
    x->amps.n = 0;
    x->amps.v = 0;
    inlet_new((t_object *)x, (t_pd *)x, gensym("list"), gensym("amps"));
    x->outlet0 = outlet_new((t_object *)x, gensym("list"));
    x->outlet1 = outlet_new((t_object *)x, gensym("list"));
    x->x_kM = 12; // The masking gradient kM is typically about 12 dB per cb
    return(x);
}

void masking_free(t_masking *x){
    freebytes(x->amps.v, x->amps.n*sizeof(float));
}

void masking_setup(void){ // (Hz, YL)->(Hz, AL)
    masking_class = class_new(gensym("masking"), (t_newmethod)masking_new,
            (t_method)masking_free, sizeof(t_masking), 0, 0);
	class_addlist(masking_class, (t_method)masking_list);
	class_addmethod(masking_class, (t_method)masking_amps, gensym("amps"), A_GIMME, 0);
	class_addmethod(masking_class, (t_method)masking_kM, gensym("kM"), A_FLOAT, 0);
}
