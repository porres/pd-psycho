// by porres, 2019, based on "Harmony: A Psychoacoustical Approach" (Parncutt, R. 1989).

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>

typedef struct floatlist{
	int n;
	t_float *v;
}floatlist;

typedef struct tonalness{
    t_object o;
    floatlist pta;
    floatlist freqs;
    floatlist amps;
    t_outlet *outlet0;
    t_outlet *outlet1;
    int ncats;
    float tolerance;
    float kM;
}tonalness;

t_class *tonalness_class;

float kT = 3;

void tonalness_kM(tonalness *x, float f){
    x->kM = f;
}

void tonalness_tolerance(tonalness *x, float f){
    x->tolerance = f;
}

float erb(float Hz){
    float kHz = Hz/1000;
    return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

static void tonalness_masking(tonalness *x, int ac){
    int i, j;
    for(i = 0; i < ac; i++){
        float kHz = x->freqs.v[i]/1000;
        float LTh = 3.64 * pow(kHz, -0.8) - 6.5 * exp(-0.6 * pow(kHz -3.3, 2))
        + 0.001 * pow(kHz, 4); // Eq.2
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
                float pthDif = fabs(erb(x->freqs.v[j]) - erb(x->freqs.v[i])); // Eq.3
                // pitch difference in critical bandwidths
                if (pthDif < 3.){ // otherwise masking is negligible
                    // Partial Masking Level (Eq.4): masking due to one masker
                    float PML = YL_j - x->kM * pthDif;
                    // The masking gradient kM is typically about 12 dB per cb
                    // Assumption: "Typical" complex tones in reg 4 have 10 audible harmonics.
                    sum += pow(10, (PML / 20.)); // add amplitudes
                }
            }
        }
        // Overall Masking Level (Eq.5)
        float ML = fmax(20.*log10(sum), 0); // due to all maskers
        // Audible Level (level above masked threshold) (Eq.6)
        float AL = fmax(YL_i-ML, 0);
        SETFLOAT(&at[i], AL);
    }
    for(i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(at+i);
}

// Tonalness

void tonalness_list(tonalness *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->freqs.v)
        freebytes(x->freqs.v, x->freqs.n*sizeof(float));
    x->freqs.v = (float *)getbytes(ac*sizeof(float));
    x->freqs.n = ac;
    int i, j, n;
    for(i = 0; i < ac; i++)
        x->freqs.v[i] = atom_getfloat(av+i);
    if(x->freqs.n != x->amps.n){
        pd_error(x, "[tonalness]: input lists don't have the same length");
        return;
    }
    tonalness_masking(x, ac);
    if(x->pta.v)
        freebytes(x->pta.v, x->pta.n*sizeof(float));
    x->pta.v = (float *)getbytes(ac*sizeof(float));
    x->pta.n = ac;
    int ncats = x->ncats;
    float PTA[ncats], CTA[ncats], MAXTA[ncats];
    for(i = 0; i < ncats; i++){
        PTA[i] = 0;
        CTA[i] = 0;
        MAXTA[i] = 0;
    }
    for(i = 0; i < ac; i++){
        float AL = x->amps.v[i];
        float thisPTA = 1.-exp(AL/-15.); // Audibility Ap(p) - Eq.7
        post("Ap(p) = %f", f);
        x->pta.v[i] =  thisPTA;
        int cat = round((ftom(x->freqs.v[i])-12) * ncats/120.); // 120.???
        if(cat >= 0 && cat < ncats)
            PTA[cat] += thisPTA;
    }
    float maxcta = 0;
    for(i = 0; i < ncats; i++){
        float sum = 0;
        float Hz_orig = mtof(12+i*120/(float)ncats);
        for(n = 1; n <= 10; n++){ // "pseudo-harmonics"
            float Hz_harm = Hz_orig * n;
            j = round((ftom(Hz_harm)-12)*ncats/120.); // 120.???
            // post("i=%d Hz_orig=%f, PTA=%f, Hz_harm=%f, PTA=%f",i,Hz_orig,PTA[i],Hz_harm,PTA[j]);
            if(j >= 0 && j < ncats)
                if(PTA[j] > 0)
                    sum += sqrt(PTA[j]/n);
        }
        CTA[i] = sum*sum / kT; // complex tone audibility (virtual pitch weight)
        if(maxcta < CTA[i])
            maxcta = CTA[i];
    }
	float sum = 0;
	for(i = 0; i < x->pta.n; i++)
        sum += (x->pta.v[i] * x->pta.v[i]);
	outlet_float(x->outlet1, maxcta/5);
	outlet_float(x->outlet0, sqrt(sum)/2);
}

void tonalness_amps(tonalness *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->amps.v)
        freebytes(x->amps.v, x->amps.n*sizeof(float));
    x->amps.v = (float *)getbytes(ac*sizeof(float));
    x->amps.n = ac;
    for(int i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(av+i);
}

void tonalness_free(tonalness *x){
    freebytes(x->pta.v, x->pta.n*sizeof(float));
    freebytes(x->freqs.v, x->freqs.n*sizeof(float));
    freebytes(x->amps.v, x->amps.n*sizeof(float));}

void *tonalness_new(t_floatarg f){
    tonalness *x = (tonalness *)pd_new(tonalness_class);
    x->pta.n = 0;
    x->pta.v = 0;
    x->freqs.n = 0;
    x->freqs.v = 0;
    x->amps.n = 0;
    x->amps.v = 0;
    inlet_new((t_object *)x, (t_pd *)x, gensym("list"), gensym("amps"));
    x->outlet0 = outlet_new((t_object *)x, gensym("list"));
    x->outlet1 = outlet_new((t_object *)x, gensym("list"));
    // 10 octaves. default is 12 pitches per octave.
    x->ncats = 10 * (f <= 0 ? 12 : f);
    x->tolerance = 0; // in cats (not used)
    x->kM = 12; // The masking gradient kM is typically about 12 dB per cb
    return(x);
}

void tonalness_setup(void){
    tonalness_class = class_new(gensym("tonalness"), (t_newmethod)tonalness_new,
            (t_method)tonalness_free, sizeof(tonalness), 0, A_DEFFLOAT, 0);
	class_addlist(tonalness_class, (t_method)tonalness_list);
    class_addmethod(tonalness_class, (t_method)tonalness_amps, gensym("amps"), A_GIMME, 0);
    class_addmethod(tonalness_class, (t_method)tonalness_tolerance,
                    gensym("tolerance"), A_FLOAT, 0);
    class_addmethod(tonalness_class, (t_method)tonalness_kM,
                    gensym("kM"), A_FLOAT, 0);
}
