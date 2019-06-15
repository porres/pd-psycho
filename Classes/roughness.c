/*
 * Copyright (C) 2019 - Alexandre Porres
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 */
// model derived from the work of Parncutt, Sethares, Vassilakis & Barlow.

#include <math.h>
#include "m_pd.h"

#define lOG10_2 log10(2)
#define EXP_1 exp(1)

static t_class *roughness_class;

typedef struct roughness{
    t_object  x_ob;
    t_outlet *x_outlet;
    float     x_curve; // 0 = parncutt / 1 = sethares
    int       x_ac;
    t_atom   *x_freqs;
    float    *x_amps;
    int      x_phonsteps;
    int      x_loudness;
    int      x_weight;
    int      x_masking;
}t_roughness;

//======================= Amplitude/Loudness functions ====================================/

/*
 static float YL(float Hz, float amp){ // ?????
 float kHz = Hz/1000;
 float LTh = 3.64 * pow(kHz,-0.8) - 6.5*exp(-0.6*pow(kHz-3.3, 2)) + 0.001*pow(kHz, 4);
 float dB = rmstodb(amp);
 return dB-LTh;
 }
 
 // ******* [dbtorms] - Converts dB to Linear Amplitude
 static float db2gain(float db){
 return pow(10, (db-100)/20);
 }
 */

// ******* [rmstodb] - Converts Linear Amplitude to dB
static float gain2db(float amp){
    if(amp > 0.00001)
        return 100+20*log10f(amp);
    return 0;
}

//***** Convert Phons and Hz to dBA (function by Clarence Barlow, data from Robson & Dadson)
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

// dbA
static float dbA(float hz, float ph){
    return dd(hz, ph) + (throfs(hz) - dd(hz, 3)) * (hz/250000);
}

//**** Convert dBA and Hz to Phons (function by Clarence Barlow, data from Robson & Dadson)

static float phon(float hz, float db, int steps){
/* NB: 1 dbA (the minimum workable value)
 gives -123.3 Ph @ 20 Hz, -30.0 Ph @ 20000 Hz */
	float ph, phbt, phtp, lphbt, lphtp;
    int i = 0;
	float mb;
	phbt = -125.95;
    phtp = 150.00;
	do{
		lphbt = phbt;
        lphtp = phtp;
		ph = (phbt+phtp)/2;
        mb = dbA(hz,ph);
        if(mb < db){
            phbt = ph;
            phtp = lphtp;
        }
        else{
            phbt = lphbt;
            phtp = ph;
        }
		i++;
	} while (i < steps);
	return ph;
}

// *************** Phons to Sones Conversion
static float ph2sn(float ph){
	if(ph >  40)
		return pow(10, lOG10_2*(ph-40)/10);
	else if(ph > 1)
		return pow(ph/40, 2.86);
	else
        return 0;
}

// *************** Masking // Used by Barlow (Masking model from Parncutt's Pitch Commonality model)
float erb (float Hz){ // equal rectangular bandwidth
	float kHz = Hz/1000;
	return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

static void do_the_masking(int ac, t_atom *freqs, float *amps){
	float a[ac];
	float kM = 18;
	int i, j;
	for(i = 0; i < ac; i++){
		float Hz = atom_getfloat(&freqs[i]);
		float kHz = Hz/1000;
		float LTh = 3.64 * pow(kHz, -0.8) -6.5 * exp(-0.6*pow(kHz-3.3, 2)) + 0.001*pow(kHz, 4);
		float dB = rmstodb(amps[i]);
		float YL_i = dB-LTh;
		float sum = 0;
		for(j = 0; j < ac; j++) if(i != j){ // masker (which doesn't mask itself)
			float YL_j = amps[j];
			float pthDif = fabs(erb(atom_getfloat(&freqs[j])) - erb(atom_getfloat(&freqs[i])));
			// pitch difference in critical bandwidths
			if(pthDif < 3.){ // (otherwise masking is negligible) // (to save computing time)
				float PML = YL_j - kM * pthDif; // masking due to one masker
				// The masking gradient kM is typically about 12 dB per cb
				// Assumption: "Typical" complex tones in reg 4 have 10 audible harmonics.
				sum += pow(10, (PML / 20.)); // add amplitudes
			}
		}
		float ML = fmax(20.*log10(sum), 0); // due to all maskers
		float AL = fmax(YL_i-ML, 0); // level above masked threshold
		a[i]=dbtorms(AL+LTh);
	}
	for(i = 0; i < ac; i++)
        amps[i] = a[i];
}

// ******** Vassilakis' weight (weight = 0)
static float vassilakis(float a1, float a2) {
    if(a1 == 0 || a2 == 0)
        return 0;
    float amin = a1 < a2 ? a1 : a2;
    float x = a1 * a2;
    float y = 2*amin/(a1+a2);
    return pow(x, 0.1) * pow(y, 3.11);
}

// ======================= Frequency to Barks ====================================/
static float barks(float hz){
    if(hz < 219.5)
        return 13.3 * atan(3*hz/4000); // Terhardt
    else{ // Traumuller
        float r = ((26.81*hz)/(1960 + hz)) - 0.53;
        if(r > 20.1)
            r = r + 0.22 * (r - 20.1);
        return r;
    }
}

// =============== Bark to Roughness (Plomp & Levelt's Curve) =====================/
static float parncutt(float bark_dif){ // Curve = 0
//    return bark_dif > 1.2 ? 0 : pow(4*EXP_1*bark_dif / exp(4*bark_dif), 2);
    return pow(4*EXP_1*bark_dif / exp(4*bark_dif), 2);
}

static float sethares(float bark_dif){ // Curve = 1
    return 5.56309 * (exp(bark_dif * -3.51) - exp(bark_dif * -5.75));
}

//************************ CALCULATE ROUGHNESS ************************/
static void roughness(t_roughness *x){
    int size = x->x_ac;
    int i, j;
    float total_roughness = 0;
    float freq_i, freq_j; // frequencies [0, 1]
    float amp_i, amp_j;   // amplitudes  [0, 1]
    float amp[size];
    for(i = 0; i < size; i++){
        freq_i = atom_getfloat(x->x_freqs + i);
        amp[i] = 0.00001 > x->x_amps[i] ? 0.00001 : x->x_amps[i];
    }
    if(x->x_masking)
        do_the_masking(size, x->x_freqs, amp);
    for(i = 0; i < size; i++){
        if(x->x_loudness == 0){ // linear amplitude (do nothing)
        } // used by Setahres and Vassilakis
        else if(x->x_loudness == 1) // Sones (used by Porres & Barlow)
            amp[i] = ph2sn(phon(freq_i, gain2db(x->x_amps[i]), x->x_phonsteps));
/*      else if(x->x_loudness == 2) // Phon
            amp[i] = phon(freq_i,gain2db(x->x_amps[i]), x->x_phonsteps);
        else if(x->x_loudness == 3) // dbtorms(dB-LTh)
            amp[i] = db2gain(YL(freq_i, gain2db(amp[i])));
        else if(x->x_loudness == 4) // dbtorms(Phon)
            amp[i] = db2gain(phon(freq_i, gain2db(x->x_amps[i]), x->x_phonsteps));
        else
            post("[roughness]: loudness mode options: 0, 1, 2, 3 & 4");*/
        else
            post("[roughness]: loudness mode options: <0> (linear) or <1> (Sones)");
    }
    for(i = 0; i < size; i++){
        for(j = i + 1; j < size; j++){
            freq_i = atom_getfloat(x->x_freqs + i);
            freq_j = atom_getfloat(x->x_freqs + j);
            amp_i = amp[i];
            amp_j = amp[j];
            float ampweight;
            if(x->x_weight == 0) // Ampweight = Vassilakis
                ampweight = vassilakis(amp_i, amp_j);
            else if(x->x_weight == 1)
                ampweight = sqrt(amp_i * amp_j); // Ampweight = Barlow
            else
                ampweight = amp_i < amp_j ? amp_i : amp_j; // Ampweight = Sethares (min)
            float bark_dif = fabs(barks(freq_i) - barks(freq_j));
            if(x->x_curve == 0) // Parncutt
                total_roughness += (parncutt(bark_dif) * ampweight);
            else // Sethares
                total_roughness += (sethares(bark_dif) * ampweight);
        }
    }
    if(total_roughness != 0 && !isnormal(total_roughness)){
        post("warning: total_roughness was %f : zeroing", total_roughness);
        startpost("x->x_amps[i] were :");
        for(i = 0; i < size; i++)
            startpost(" %f", x->x_amps[i]);
        post("");
        startpost("amp[i] were :");
        for (i=0; i<size; i++) startpost(" %f", amp[i]);
        post("");
        total_roughness = 0;
    }
    outlet_float(x->x_outlet, total_roughness);
}

// =============== METHODS =====================/
// Frequency input
static void roughness_freqs(t_roughness *x, t_symbol *s, int ac, t_atom *av){
    // LATER accelerqte x->x_freqs by making it a float *.
    // LATER see if resizebytes gets any faster
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->x_freqs)
        freebytes(x->x_freqs, x->x_ac*sizeof(t_atom));
    x->x_freqs = (t_atom *)getbytes(ac * sizeof(t_atom));
    for(int i = 0; i < ac; i++)
        x->x_freqs[i] = av[i];
    x->x_ac = ac;
    if(x->x_amps)
        roughness(x);
}

// Amplitude input
static void roughness_amps(t_roughness *x, t_symbol *s, int ac, t_atom *av){
    // LATER see if resizebytes gets any faster
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->x_amps)
        freebytes(x->x_amps,x->x_ac*sizeof(float));
    x->x_amps = (float *)getbytes(ac * sizeof(float));
	x->x_ac = ac;
    for(int i = 0; i < ac;i++)
		x->x_amps[i] = atom_getfloat(av+i);
}

static void roughness_phonsteps(t_roughness *x, t_float f){ // undocumented, for tests only
    x->x_phonsteps = f;
}

static void roughness_loudness(t_roughness *x, t_float f){
        x->x_loudness = (int)f != 0;
}

static void roughness_weight(t_roughness *x, t_float f){
    int i = (int)f;
    x->x_weight = i < 0 ? 0 : i > 2 ? 2 : i;
}

static void roughness_curve(t_roughness *x, t_float f){
    x->x_curve = f != 0;
}

static void roughness_masking(t_roughness *x, t_float f){
	x->x_masking = f != 0;
}

// Models
static void roughness_barlow(t_roughness *x){ // Clarence Barlow
    x->x_loudness = 1;
    x->x_masking = 1;
    x->x_curve = 0;
    x->x_weight = 1;
}

static void roughness_porres(t_roughness *x){ // Alexandre Porres
    x->x_loudness = 1;
    x->x_masking = 1;
    x->x_curve = 0;
    x->x_weight = 0;
}

static void roughness_vass(t_roughness *x){ // Pantelis Vassilakis
    x->x_loudness = 0;
    x->x_masking = 0;
    x->x_curve = 1;
    x->x_weight = 0;
}

static void roughness_seth(t_roughness *x){ // William Sethares
    x->x_loudness = 0;
    x->x_masking = 0;
    x->x_curve = 1;
    x->x_weight = 2;
}

// Free
static void roughness_free (t_roughness *x){
    if(x->x_freqs)
        freebytes(x->x_freqs, x->x_ac*sizeof(t_atom));
    if(x->x_amps)
        freebytes(x->x_amps, x->x_ac*sizeof(float));
}

// New
static void *roughness_new(void){
    t_roughness *x = (t_roughness *)pd_new(roughness_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("list"), gensym("amps"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->x_phonsteps = 20; // hardcoded, can be changed, but it's undocumented
    roughness_porres(x); // default
    return(void *)x;
}

void roughness_setup(void){
    roughness_class = class_new(gensym("roughness"), (t_newmethod)roughness_new, (t_method)roughness_free, sizeof(t_roughness), 0, 0);
	class_addlist(roughness_class, roughness_freqs);
	class_addmethod(roughness_class, (t_method)roughness_amps,
                    gensym("amps"), A_GIMME, 0);
    class_addmethod(roughness_class, (t_method)roughness_phonsteps,
                    gensym("phonsteps"), A_FLOAT, 0); // undocumented
    class_addmethod(roughness_class, (t_method)roughness_loudness,
                    gensym("loudness"), A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_weight,
                    gensym("weight"), A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_curve,
                    gensym("curve"), A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_masking,
                    gensym("masking"), A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_seth,
                    gensym("sethares"), 0);
    class_addmethod(roughness_class, (t_method)roughness_vass,
                    gensym("vassilakis"), 0);
    class_addmethod(roughness_class, (t_method)roughness_barlow,
                    gensym("barlow"), 0);
    class_addmethod(roughness_class, (t_method)roughness_porres,
                    gensym("porres"), 0);
}
