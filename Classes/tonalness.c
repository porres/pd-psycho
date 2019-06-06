// by porres, 2019, based on "Harmony: A Psychoacoustical Approach" (Parncutt, R. 1989).

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>

typedef struct floatlist{
	int         n;
	t_float    *v;
}t_floatlist;

typedef struct tonalness{
    t_object    x_obj;
    t_floatlist pta;
    t_floatlist amps;
    t_outlet   *outlet0;
    t_outlet   *outlet1;
    int         x_ncats;
    float       x_tolerance;
    float       x_kM;
}t_tonalness;

t_class *tonalness_class;

float kT = 3;

void tonalness_kM(t_tonalness *x, float f){
    x->x_kM = f;
}

void tonalness_tolerance(t_tonalness *x, float f){
    x->x_tolerance = f;
}

float erb(float Hz){
    float kHz = Hz/1000;
    return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

// Tonalness

void tonalness_list(t_tonalness *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(ac != x->amps.n){
        pd_error(x, "[tonalness]: input lists don't have the same length");
        return;
    }
    int ncats = x->x_ncats;
    float PTA[ncats], CTA[ncats];
    int i, j, n;
    float sum_pta = 0, maxcta = 0, sum = 0;
    for(i = 0; i < ncats; i++){
        PTA[i] = 0;
        CTA[i] = 0;
    }
    for(int i = 0; i < ac; i++){ // YL
        float kHz = atom_getfloat(av+i)/1000;
        float LTh = 3.64 * pow(kHz, -0.8) - 6.5 * exp(-0.6 * pow(kHz -3.3, 2))
        + 0.001 * pow(kHz, 4); // Eq.2
        float dB = rmstodb(x->amps.v[i]);
        x->amps.v[i] = dB-LTh; // YL
    }
    for(i = 0; i < ac; i++){ // PTA
        float YL_i = x->amps.v[i];
        sum = 0;
        for(j = 0; j < ac; j++){ // MASKING
            if(i != j){
                float YL_j = x->amps.v[j];
                float dif = fabs(erb(atom_getfloat(av+j)) - erb(atom_getfloat(av+i)));
                // pitch difference in critical bandwidths
                if(dif < 3.){ // otherwise masking is negligible
                    float PML = YL_j - x->x_kM * dif; // (Eq.4) Partial Masking Level
                    // The masking gradient kM is typically about 12 dB per cb
                    // Assumption: "Typical" complex tones in reg 4 have 10 audible harmonics.
                    sum += pow(10, PML/20.);
                }
            }
        }
        float ML = fmax(20.*log10(sum), 0); // (Eq.5) Masking Level
        float AL = fmax(YL_i-ML, 0); // (Eq.6) Audible Level (above masked threshold)
        float thisPTA = 1.-exp(-AL/15.); // (Eq.7) Audibility Ap(P)
        sum_pta += (thisPTA * thisPTA);;
        int cat = round((ftom(atom_getfloat(av+i))-12)); // -12???
        if(cat >= 0 && cat < ncats) // why?
            PTA[cat] += thisPTA;
    }
    for(i = 0; i < ncats; i++){ // CTA
        sum = 0;
        float Hz_orig = mtof(12+i); // + 12???
        for(n = 1; n <= 10; n++){ // "pseudo-harmonics"
            float Hz_harm = Hz_orig * n;
            j = round((ftom(Hz_harm)-12)); // -12???
            if(j >= 0 && j < ncats && PTA[j] > 0)
                sum += sqrt(PTA[j]/n);
        }
        CTA[i] = sum*sum / kT; // complex tone audibility (virtual pitch weight)
        if(maxcta < CTA[i])
            maxcta = CTA[i];
    }
	outlet_float(x->outlet1, maxcta*0.2);
	outlet_float(x->outlet0, sqrt(sum_pta)*0.5);
}

void tonalness_amps(t_tonalness *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->amps.v)
        freebytes(x->amps.v, x->amps.n*sizeof(float));
    x->amps.v = (float *)getbytes(ac*sizeof(float));
    x->amps.n = ac;
    for(int i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(av+i);
}

void tonalness_free(t_tonalness *x){
    freebytes(x->pta.v, x->pta.n*sizeof(float));
    freebytes(x->amps.v, x->amps.n*sizeof(float));}

void *tonalness_new(t_floatarg f){
    t_tonalness *x = (t_tonalness *)pd_new(tonalness_class);
    x->amps.n = 0;
    x->amps.v = 0;
    inlet_new((t_object *)x, (t_pd *)x, gensym("list"), gensym("amps"));
    x->outlet0 = outlet_new((t_object *)x, gensym("list"));
    x->outlet1 = outlet_new((t_object *)x, gensym("list"));
    // 10 octaves. default is 12 pitches per octave.
    x->x_ncats = 10 * (f <= 0 ? 12 : f);
    x->x_tolerance = 0; // in cats (not used)
    x->x_kM = 12; // The masking gradient kM is typically about 12 dB per cb
    return(x);
}

void tonalness_setup(void){
    tonalness_class = class_new(gensym("tonalness"), (t_newmethod)tonalness_new,
            (t_method)tonalness_free, sizeof(t_tonalness), 0, A_DEFFLOAT, 0);
	class_addlist(tonalness_class, (t_method)tonalness_list);
    class_addmethod(tonalness_class, (t_method)tonalness_amps, gensym("amps"), A_GIMME, 0);
    class_addmethod(tonalness_class, (t_method)tonalness_tolerance,
                    gensym("tolerance"), A_FLOAT, 0);
    class_addmethod(tonalness_class, (t_method)tonalness_kM,
                    gensym("kM"), A_FLOAT, 0);
}
