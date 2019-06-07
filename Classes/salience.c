/*
 By porres, 2010-2019, based on:
 Parncutt, Richard & Strasburger, Hans. (1994).
 Applying psychoacoustics in composition: Harmonic progressions of non-harmonic sonorities.
 Perspectives of New Music. 32. 88-129. 10.2307/833600.
 */

#include "m_pd.h"
#include <math.h>

typedef struct floatlist{
	int n;
	t_float *v;
}floatlist;

// Pitch Salience Profile -------------------------------------------------------------

typedef struct salience{
	t_object    o;
    floatlist   freqs;
    floatlist   amps;
	t_outlet   *outlet0;
	t_outlet   *outlet1;
    int         ncats;
    float       tolerance;
    float       kM;
}salience;

static t_class *salience_class;

float kT = 3;

// Masking

float erb(float Hz){
    float kHz = Hz/1000;
    return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

static void salience_masking(salience *x, int ac){
    int i, j;
    for(i = 0; i < ac; i++){
        float kHz = x->freqs.v[i]/1000;
        float LTh = 3.64 * pow(kHz, -0.8) - 6.5 * exp(-0.6 * pow(kHz -3.3, 2))
        + 0.001 * pow(kHz, 4);
        float dB = rmstodb(x->amps.v[i]);
        x->amps.v[i] = dB-LTh;
    }
    t_atom at[ac];
    for(i = 0; i < ac; i++){
        float YL_i = x->amps.v[i];
        float sum = 0;
        for(j = 0; j < ac; j++){
            if(i != j){ // masker (which doesn't mask itself)
                float YL_j = x->amps.v[j];
                float pthDif = fabs(erb(x->freqs.v[j]) - erb(x->freqs.v[i]));
                // pitch difference in critical bandwidths
                if (pthDif < 3.){ // otherwise masking is negligible
                    float PML = YL_j - x->kM * pthDif; // masking due to one masker
                    // The masking gradient kM is typically about 12 dB per cb
                    // Assumption: "Typical" complex tones in reg 4 have 10 audible harmonics.
                    sum += pow(10, (PML / 20.)); // add amplitudes
                }
            }
        }
        float ML = fmax(20.*log10(sum), 0); // due to all maskers
        float AL = fmax(YL_i-ML, 0); // AL: level above masked threshold
        SETFLOAT(&at[i], AL);
    }
    for(i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(at+i);
}

// Methods

void salience_tolerance(salience *x, float f){
    x->tolerance = f;
}

void salience_kM(salience *x, float f){
    x->kM = f;
}

void salience_list(salience *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    int i, j, n;
    float max_ta = 0;
    float total_max_ta = 0;
    if(x->freqs.v)
        freebytes(x->freqs.v, x->freqs.n*sizeof(float));
    x->freqs.v = (float *)getbytes(ac*sizeof(float));
    x->freqs.n = ac;
    for(i = 0; i < ac; i++)
        x->freqs.v[i] = atom_getfloat(av+i);
    if(x->freqs.n != x->amps.n){
        pd_error(x, "[salience]: input lists don't have the same length");
        return;
    }
    salience_masking(x, ac);
// weigth
    int ncats = x->ncats;
    float PTA[ncats], CTA[ncats], MAXTA[ncats];
    for(i = 0; i < ncats; i++){
        PTA[i] = 0;
        CTA[i] = 0;
        MAXTA[i] = 0;
    }
    for(i = 0; i < ac; i++){
        float AL = x->amps.v[i];
        float thisPTA = 1.-exp(AL/-15.);
        int cat = round((ftom(x->freqs.v[i])-12) * ncats/120.); // 120.???
        if(cat >= 0 && cat < ncats)
            PTA[cat] += thisPTA;
    }
    for(i = 0; i < ncats; i++){
        float sum = 0;
        float Hz_orig = mtof(12+i*120/(float)ncats);
        for(n = 1; n <= 10; n++){ // "pseudo-harmonics"
            float Hz_harm = Hz_orig * n;
            j = round((ftom(Hz_harm)-12)*ncats/120.); // 120.???
            if(j >= 0 && j < ncats)
                if(PTA[j] > 0)
                    sum += sqrt(PTA[j]/n);
        }
        CTA[i] = sum*sum / kT; // complex tone audibility (virtual pitch weight)
        MAXTA[i] = fmax(PTA[i], CTA[i]);
        if(max_ta < MAXTA[i])
            max_ta = MAXTA[i];
        total_max_ta += MAXTA[i];
    }
// Calculate multiplicity (number of simultaneously noticed tones or pitches)
	float mulX = total_max_ta / max_ta;
	float mul = pow(mulX, 0.5); // kS is simultaneity perception parameter, typically 0.5
	t_atom salience[ncats];
	// M - Calculate tone (pitch) salience arrays
	for(i = 0; i < ncats; i++)
        SETFLOAT(salience+i, (MAXTA[i]/max_ta) * (mul/mulX));
	outlet_float(x->outlet1, mul);
	outlet_list(x->outlet0, &s_list, ncats, salience);
}

void salience_amps(salience *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->amps.v)
        freebytes(x->amps.v, x->amps.n*sizeof(float));
    x->amps.v = (float *)getbytes(ac*sizeof(float));
    x->amps.n = ac;
    for(int i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(av+i);
}

void *salience_new(t_floatarg f){
    salience *x = (salience *)pd_new(salience_class);
    x->freqs.n = 0;
    x->freqs.v = 0;
    x->amps.n = 0;
    x->amps.v = 0;
    inlet_new((t_object *)x, (t_pd *)x, gensym("list"), gensym("amps"));
    x->outlet0 = outlet_new((t_object *)x, gensym("list"));
    x->outlet1 = outlet_new((t_object *)x, gensym("float"));
// 10 octaves. default is 12 pitches per octave.
    x->ncats = 10 * (f <= 0 ? 12 : f);
    x->tolerance = 0; // in cats (not used)
    x->kM = 18;
    return(x);
}

void salience_setup(void){ // (MAXTA[cat]) -> (SalTA, mul)
    salience_class = class_new(gensym("salience"), (t_newmethod)salience_new,
            0, sizeof(salience), 0, A_DEFFLOAT, 0);
	class_addlist(salience_class, (t_method)salience_list);
    class_addmethod(salience_class, (t_method)salience_amps,
                    gensym("amps"), A_GIMME, 0);
    class_addmethod(salience_class, (t_method)salience_tolerance,
                    gensym("tolerance"), A_FLOAT, 0);
    class_addmethod(salience_class, (t_method)salience_kM,
                    gensym("kM"), A_FLOAT, 0);
}
