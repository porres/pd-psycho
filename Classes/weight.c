// by porres, 2019, based on "Harmony: A Psychoacoustical Approach" (Parncutt, R. 1989).

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>

typedef struct floatlist{
	int         n;
	t_float    *v;
}t_floatlist;

typedef struct weight{
    t_object    x_obj;
    t_floatlist amps;
    t_outlet   *outlet0;
    t_outlet   *outlet1;
    int         x_ncats;
    float       x_tolerance;
    float       x_kM;
    t_outlet   *outlet2;
}t_weight;

t_class *weight_class;

float kT = 3;

// Masking

void weight_kM(t_weight *x, float f){
    x->x_kM = f;
}

float erb(float Hz){
    float kHz = Hz/1000;
    return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

static void weight_masking(t_weight *x, int ac, t_atom *av){
    int i, j;
    for(i = 0; i < ac; i++){
        float kHz = atom_getfloat(av+i)/1000;
        float LTh = 3.64 * pow(kHz, -0.8) - 6.5 * exp(-0.6 * pow(kHz -3.3, 2))
        + 0.001 * pow(kHz, 4); // Eq.2
        float dB = rmstodb(x->amps.v[i]);
        x->amps.v[i] = dB-LTh; // YL
    }
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
        x->amps.v[i] = fmax(YL_i-ML, 0);
    }
}

void weight_list(t_weight *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(ac != x->amps.n){
        pd_error(x, "[weight]: input lists don't have the same length");
        return;
    }
    weight_masking(x, ac, av);
	t_atom PTAout[ac];
	int ncats = x->x_ncats;
	float PTA[ncats], CTA[ncats], MAXTA[ncats];
	int i, j, n;
	for(i = 0; i < ncats; i++){
		PTA[i] = 0;
		CTA[i] = 0;
		MAXTA[i] = 0;
	}
	for(i = 0; i < ac; i++){
		float AL = x->amps.v[i];
//        post("AL = %f", AL);
		float thisPTA = 1.-exp(AL/-15.);
//        post("thisPTA = %f", thisPTA);
		SETFLOAT(&PTAout[i], thisPTA);
		int cat = round((ftom(atom_getfloat(av+i))-12)); // -12???
		if(cat >= 0 && cat < ncats) // why?
            PTA[cat] += thisPTA;
	}
	for(i = 0; i < ncats; i++){
		float sum = 0;
		float Hz_orig = mtof(12+i*120/(float)ncats);
		for(n = 1; n <= 10; n++){ // "pseudo-harmonics"
			float Hz_harm = Hz_orig * n;
			j = round((ftom(Hz_harm)-12)); // -12???
            if(j >= 0 && j < ncats && PTA[j] > 0)
                    sum += sqrt(PTA[j]/n);
		}
		CTA[i] = sum*sum / kT; // complex tone audibility (virtual pitch weight)
		MAXTA[i] = fmax(PTA[i], CTA[i]);
	}
    t_atom at[ncats];
    for(i = 0; i < ncats; i++)
        SETFLOAT(at+i, MAXTA[i]);
    outlet_list(x->outlet2, &s_list, ncats, at);
    for(i = 0; i < ncats; i++)
        SETFLOAT(at+i, CTA[i]);
    outlet_list(x->outlet1, &s_list, ncats, at);
	outlet_list(x->outlet0, &s_list, ac, PTAout);
}

void weight_amps(t_weight *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->amps.v)
        freebytes(x->amps.v, x->amps.n*sizeof(float));
    x->amps.v = (float *)getbytes(ac*sizeof(float));
    x->amps.n = ac;
    for(int i = 0; i < ac; i++)
        x->amps.v[i] = atom_getfloat(av+i);
}

void weight_tolerance(t_weight *x, float f){
    x->x_tolerance = f;
}

void weight_free(t_weight *x){
    freebytes(x->amps.v, x->amps.n*sizeof(float));
}

void *weight_new(t_floatarg f){
    t_weight *x = (t_weight *)pd_new(weight_class);
    x->amps.n = 0;
    x->amps.v = 0;
    inlet_new((t_object *)x, (t_pd *)x, gensym("list"), gensym("amps"));
    x->outlet0 = outlet_new((t_object *)x, gensym("list"));
    x->outlet1 = outlet_new((t_object *)x, gensym("list"));
    x->outlet2 = outlet_new((t_object *)x, gensym("list"));
    // 10 octaves. default is 12 pitches per octave.
    x->x_ncats = 10 * (f <= 0 ? 12 : f);
    x->x_tolerance = 0; // in cats (not used)
    x->x_kM = 12;
    return(x);
}

void weight_setup(void){ // (Hz,AL) -> (Hz, PTA, CTA[cat], MAXTA[cat])
    weight_class = class_new(gensym("weight"), (t_newmethod)weight_new,
            (t_method)weight_free, sizeof(t_weight), 0, A_DEFFLOAT, 0);
	class_addlist(weight_class, (t_method)weight_list);
	class_addmethod(weight_class, (t_method)weight_amps,
                    gensym("amps"), A_GIMME, 0);
	class_addmethod(weight_class, (t_method)weight_tolerance,
                    gensym("tolerance"), A_FLOAT, 0);
    class_addmethod(weight_class, (t_method)weight_kM,
                    gensym("kM"), A_FLOAT, 0);
}
