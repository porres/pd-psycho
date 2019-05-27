/* Psycho library
 *
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
 *
 */

#include <string.h>
#include <math.h>
#include "m_pd.h"

static t_class *roughness_class;

typedef enum bool {false, true} bool;

typedef struct roughness{
    t_object x_ob;
    t_outlet *x_outlet;
	//config
	float curve; //0 = parncutt / 1 = sethares (??????????????)
	int argc;
    t_atom *freqs;
	float *amps;
	int phonsteps;
	int loudness;
	int mean; // (!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
	int masking;
}t_roughness;


//________FUNCTIONS__________________________________________________________________/

// GET RID OF THESE!!!!

// ************************ min
static float min(float a, float b){
  return a < b ? a : b;
}

// ************************ max
static float max(float a, float b){
  return a > b ? a : b;
}

// ************************ Square

static float sqr(float v){
    return(v*v);
}

//!!!!! PSYCHOACOUSTICAL FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/

//=============================== Amplitude & Loudness ====================================/

// ******* [rmstodb] - Converts Linear Amplitude to dB(A)
static float gain2db(float amp) {
    if (amp > 0.00001)
        return 100+20*log10f(amp);
    return 0;
}

// ******* [dbtorms] - Converts dB(A) to Linear Amplitude
static float db2gain(float db){
    return pow(10, (db-100)/20);
}

// ******* Convert dBA and Hz to Phons (function by Clarence Barlow, data from Robson & Dadson)

static float dbA(float hz, float ph);

static float phon(float hz, float db, int phonsteps)
/* NB: 1 dbA (the minimum workable value)
 gives -123.3 Ph @ 20 Hz, -30.0 Ph @ 20000 Hz */
{
	float ph,phbt,phtp,lphbt,lphtp; int i;
	float mb;
	phbt = -125.95; phtp = 150.00; i = 0;
	do{
		lphbt = phbt; lphtp = phtp;
		ph = (phbt+phtp)/2; mb = dbA(hz,ph);
		if (mb<db) { phbt = ph; phtp = lphtp; }
		else       { phbt = lphbt; phtp = ph; }
		i++;
	} while (i < phonsteps);
	return ph;
}

//**************** Convert Phons and Hz to dBA (function by Clarence Barlow, data from Robson & Dadson)

static float throfs(float hz){
    float khz, thr;
	khz = hz/1000;
	thr = 3.64*exp(-0.8*log(khz))+0.001*sqr(sqr(khz));
	if(khz <= 15.41)
        thr = thr-6.5*exp(-0.6*sqr(khz-3.3));
	return(thr);
}

float e(float a){
    if (fabs(a) < 80)
        return(exp(a));
    else return(0);
}

float bell(float a, float b){
    if (b != 0)
        return(e(-sqr(a)/(b*16)));
    else return(0);
}

float tanhyp(float x){
    float u = e(x);
    float v = e(-x);
    if ((u+v) != 0)
        return((u-v)/(u+v));
    else
        return(0);
}

float hypt(float a, float b){
    if(b != 0)
        return (1-tanhyp(a/b));
    else
        return (0);
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
    return(rbf);
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
	return(dd (hz, ph) + (throfs(hz) - dd(hz,3)) * (hz/250000));
}

 // *************** Phons to Sones Conversion
 //    0.30103 = log10(2)
static float ph2sn(float ph){
	if (ph >  40)
		return pow(10, (0.30103*(ph-40))/10);
	else if (ph > 1)
		return pow(ph/40,2.86);
	else return 0;
}

/*
 ************************ Sones to Phons Convertion

// unused
 
static float sn2ph(float sn) {
    return 40 + 33.22*log10f(sn); // fdgdf
} */

//==================================== Amplitude Weight Functions ====================================/

/*
 ************************ Vassilakis (Ampweight = 0)
 */

static float vassilakis(float a1, float a2) {
	if (a1==0 || a2==0) return 0; //???
	float amin = min(a1,a2);
	float x = a1 * a2;
	float y = 2*amin/(a1+a2);
	return (pow(x,0.1) * 0.5 * pow(y,3.11)*2);
}

//==================================== Frequency to Barks ====================================/

//***************** Terhardt
static float terh(float h){
	return 13.3*atan(3*h/4000);
}

// ******************** Traumuller
static float traun(float h) {
	float r = ((26.81*h)/(1960 + h)) - 0.53;
	if (r > 20.1) {
		r = r + 0.22*(r - 20.1);
	}
	return r;
}

// ****************** Hz to Bk
static float barks(float hz){
    if (hz < 219.5)
        return terh(hz);
    else
        return traun(hz);
}

 //====================== Bark to Roughness (Plomp & Levelt) ==========================/


// ************************ Parncutt (Curve = 0)
static float parncutt(float freq1, float freq2){
    float r = barks(freq1) - barks(freq2);
    r = fabs(r);
//    if (r > 1.2) {
//        r = 0;
//    } else {
        r = pow((4*exp(1)*r)/(exp(4*r)), 2);
//    }
    return r;
}

// ************************ Sethares (Curve = 1)
static float sethares(float freq1, float freq2){
    float r = barks(freq1) - barks(freq2);
    r = fabs(r);
    r = 5.56309 * (exp(r * -3.51) - exp(r * -5.75));
    return r;
}

//******* ************************ ************************ ************************ ************************ ************************/
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ROUGHNESS FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/
//*******  ************************ ************************ ************************ ************************ ***********************/

float erb (float Hz) {
	float kHz = Hz/1000;
	return 11.17 * log((kHz + 0.312) / (kHz + 14.675)) + 43.0;
}

static float YL (float Hz, float amp) {
	float kHz = Hz/1000;
	float LTh = 3.64*pow(kHz,-0.8)-6.5*exp(-0.6*pow(kHz-3.3,2))+0.001*pow(kHz,4);
	float dB = rmstodb(amp);
	return dB-LTh;
}

static void do_the_masking (int argc, t_atom *freqs, float *amps) {
	float a[argc];
	float kM = 18;
	int i,j;
	for (i=0; i<argc; i++) {
		float Hz = atom_getfloat(&freqs[i]);
		float kHz = Hz/1000;
		float LTh = 3.64*pow(kHz,-0.8)-6.5*exp(-0.6*pow(kHz-3.3,2))+0.001*pow(kHz,4);
		float dB = rmstodb(amps[i]);
		float YL_i = dB-LTh;
		float sum = 0;
		for (j=0; j<argc; j++) if (i!=j) { // masker (which doesn't mask itself)
			float YL_j = amps[j];
			float pthDif = fabs(erb(atom_getfloat(&freqs[j])) - erb(atom_getfloat(&freqs[i])));
			// pitch difference in critical bandwidths
			if (pthDif < 3.) { // (otherwise masking is negligible) // (to save computing time)
				float PML = YL_j - kM * pthDif; // masking due to one masker
				// The masking gradient kM is typically about 12 dB per cb
				// Assumption: "Typical" complex tones in reg 4 have 10 audible harmonics.
				sum += pow(10, (PML / 20.)); // add amplitudes
			}
		}
		float ML = fmax(20.*log10(sum),0); // due to all maskers
		float AL = fmax(YL_i-ML,0); // level above masked threshold
		a[i]=dbtorms(AL+LTh);
	}
	for (i=0; i<argc; i++) amps[i]=a[i];
}

static void roughness(t_roughness *x) {
	int size = x->argc;
	int i, j;
	float froughness = 0;
	float r;
	float freq_i, freq_j; //frequencies [0,1]
	float amp_i, amp_j;   //amplitudes  [0,1]
	float amp[size];
	
	//post("using phonsteps=%d loudness=%d",x->phonsteps,x->loudness);
	//startpost("x->amps = [");
	
	for (i=0; i<size; i++) {
		//startpost("%f%c",x->amps[i],i==size-1?']':',');
		freq_i = atom_getfloat(x->freqs + i);
		amp[i] = max(0.00001,x->amps[i]);
	}
	if (x->masking) do_the_masking(size,x->freqs,amp);
	for (i=0; i<size; i++) {
		// loudness 0 : linear amplitude
		// loudness 1 : dbtorms(dB-LTh)
		// loudness 2 : dbtorms(Phon)
		// loudness 3 : Phon
		// loudness 4 : Sones
		if      (x->loudness==1) amp[i] = db2gain(  YL(freq_i,gain2db(amp[i])));
		else if (x->loudness==2) amp[i] = db2gain(phon(freq_i,gain2db(x->amps[i]),x->phonsteps));
		else if (x->loudness==3) amp[i] =         phon(freq_i,gain2db(x->amps[i]),x->phonsteps) ;
		else if (x->loudness==4) amp[i] =   ph2sn(phon(freq_i,gain2db(x->amps[i]),x->phonsteps));
		else post("only 0,1,2,3,4");
	}

    //post("");
	
	for (i = 0; i < size; i++) {
		for (j = i + 1; j < size; j++) {
            freq_i = atom_getfloat(x->freqs + i);
            freq_j = atom_getfloat(x->freqs + j);
			amp_i = amp[i];
			amp_j = amp[j];
			float ampweight;
			if (x->mean == 0)                                  //Ampweight = Vassilakis
				ampweight = vassilakis(amp_i,amp_j);           
			else
				ampweight = sqrt(amp_i*amp_j);		      //Ampweight = Barlow	
			
			if (x->curve == 0) {                               //Parncutt for Plomp & Levelt
        		r = parncutt(freq_i, freq_j) * ampweight;
				froughness += r;
			} else if (x->curve == 1) {                        //Sethares for Plomp & Levelt
				r = sethares(freq_i, freq_j) * ampweight;
				froughness += r;
			}
		}
	}
	if (froughness!=0 && !isnormal(froughness)) {
		post("warning: froughness was %f : zeroing",froughness);
		startpost("x->amps[i] were :");
		for (i=0; i<size; i++) startpost(" %f",x->amps[i]);
		post("");
		startpost("amp[i] were :");
		for (i=0; i<size; i++) startpost(" %f",amp[i]);
		post("");
	    froughness=0;
	}
	outlet_float(x->x_outlet, froughness);
}

//______END OF FUNCTIONS__________________________________________________________________/


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OBJECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//************************ ROUGHNESS OBJECT ************************/

// Roughness modelling is Based on the work of Parncutt, Sethares, Vassilakis & Clarence Barlow.

/**
 * Hot inlet
 * Allocate frequency values and execute roughness method.
 */
static void roughness_freqs(t_roughness *x, t_symbol *s, int argc, t_atom *argv)
{
    // LATER accelerqte x->freqs by making it a float *.
    // LATER see if resizebytes gets any faster
    t_symbol *dummy = s;
    dummy = NULL;
    if (x->freqs) freebytes(x->freqs,x->argc*sizeof(t_atom));
    x->freqs = (t_atom *)getbytes(argc * sizeof(t_atom));
    int i;
    for (i=0; i<argc; i++) x->freqs[i] = argv[i];
    x->argc = argc;
    if (x->amps) {
        roughness(x);
    }
}

/**
 * Cold inlet
 * Allocate amplitude values.
 */
static void roughness_amps(t_roughness *x, t_symbol *s, int argc, t_atom *argv)
{
    // LATER see if resizebytes gets any faster
    t_symbol *dummy = s;
    dummy = NULL;
    if (x->amps ) freebytes(x->amps,x->argc*sizeof(float));
    x->amps = (float *)getbytes(argc * sizeof(float));
	x->argc = argc;
	int i;
    for (i = 0; i < argc;i++) {
		x->amps[i] = atom_getfloat(argv+i);
	}
}

static void roughness_phonsteps (t_roughness *x, t_float f) { x->phonsteps = f; }
static void roughness_loudness   (t_roughness *x, t_float f) {
	int n = (int)f;
	if (n<0 || n>4) {
		pd_error(x,"loudness is supposed to be one of 0,1,2,3,4 but got %f",f);
	} else x->loudness = f;
}

static void roughness_mean (t_roughness *x, t_float f) {
	if (f!=0 && f!=1) {
		pd_error(x,"weight is supposed to be 0 or 1, but got %f",f);
	} else x->mean = f;
}

static void roughness_curve (t_roughness *x, t_float f) {
	if (f!=0 && f!=1) {
		pd_error(x,"curve is supposed to be 0 or 1, but got %f",f);
	} else x->curve = f;
}

static void roughness_masking (t_roughness *x, t_float f) {
	if (f!=0 && f!=1) {
		pd_error(x,"masking is supposed to be 0 or 1, but got %f",f);
	} else x->masking = f;
}

/*
 * New method
 */
static void *roughness_new(t_symbol *s, int argcount, t_atom *argvec)
{
    t_roughness *x = (t_roughness *)pd_new(roughness_class);
    t_symbol *dummy = s;
    dummy = NULL;
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("list"), gensym("amps"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
	//these are the valid choices
    x->curve = 0;
    x->phonsteps = 20;
	x->loudness = 1;
	x->mean = 0;
	x->masking = 0;
    if (argcount == 0 || (argvec[0].a_w.w_float != 0 && argvec[0].a_w.w_float != 1 && argvec[0].a_w.w_float != 2)) {
		post("Missing valid argument:  roughness [0|1]. Using default: '0 = parncutt'.");
	} else {
        x->curve = argvec[0].a_w.w_float;
    }
    if (argcount > 1) {
        post("warning: too many arguments");
    }
    return (void *)x;
}

static void roughness_free (t_roughness *x) {
    if (x->freqs) freebytes(x->freqs,x->argc*sizeof(t_atom));
    if (x->amps ) freebytes(x->amps ,x->argc*sizeof(float));
}


//************************ FLUNSON OBJECT ************************/

t_class *flunson_class;

typedef struct flunson {
    t_object x_ob;
    t_outlet *x_outlet;
	int argc;
	float *amps;
} t_flunson;

void *flunson_new(void){
    t_flunson *x = (t_flunson *)pd_new(flunson_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("list"), gensym("amps"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    return (void *)x;
}

void flunson_freqs(t_flunson *x, t_symbol *s, int argc, t_atom *argv)
{
    t_symbol *dummy = s;
    dummy = NULL;
    if (argc > x->argc) argc = x->argc; // too many args.
    t_atom outs[argc];
    int i;
//  for (i=0; i<argc; i++) SETFLOAT(outs+i,        phon(atom_getfloat(argv+i),        x->amps[i] ) );
    for (i=0; i<argc; i++) SETFLOAT(outs+i,db2gain(phon(atom_getfloat(argv+i),gain2db(x->amps[i]),20)));
//  for (i=0; i<argc; i++) SETFLOAT(outs+i,db2gain( dbA(atom_getfloat(argv+i),gain2db(x->amps[i]))));
    outlet_list(x->x_outlet,&s_list,argc,outs);
}

void flunson_amps(t_flunson *x, t_symbol *s, int argc, t_atom *argv)
{
    t_symbol *dummy = s;
    dummy = NULL;
    // LATER see if resizebytes gets any faster
    if (x->amps ) freebytes(x->amps,x->argc*sizeof(float));
    x->amps = (float *)getbytes(argc * sizeof(float));
	x->argc = argc;
	int i;
    for (i = 0; i < argc;i++) {
		x->amps[i] = atom_getfloat(argv+i);
	}
}

void flunson_free (t_flunson *x) {
    if (x->amps ) freebytes(x->amps ,x->argc*sizeof(float));
}

//************************ Phon / dBA / iso226 / iso226b ************************/

static t_class *phon_class;
static t_class *dbA_class;
static t_class *iso226_class; /* not a t_func */
static t_class *iso226b_class;
typedef struct func {
    t_object x_ob;
    t_outlet *x_outlet;
	float right;
} t_func;

static void *phon_new(t_float right)
{
    t_func *x = (t_func *)pd_new(phon_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("float"), gensym("right"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->right = right;
    return (void *)x;
}

static void *dbA_new(t_float right)
{
    t_func *x = (t_func *)pd_new(dbA_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("float"), gensym("right"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->right = right;
    return (void *)x;
}

static void *iso226b_new(t_float right)
{
    t_func *x = (t_func *)pd_new(iso226b_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("float"), gensym("right"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->right = right;
    return (void *)x;
}
static void function_right (t_func *x, float right) {x->right = right;}
static void   phon_left (t_func *x, float left) {outlet_float(x->x_outlet,  phon(left,x->right,20));}
static void    dbA_left (t_func *x, float left) {outlet_float(x->x_outlet,   dbA(left,x->right));}

typedef struct iso226 {
        t_object x_ob;
        t_outlet *x_outlet;
} t_iso226;

static void *iso226_new(void){
    t_func *x = (t_func *)pd_new(iso226_class);
    x->x_outlet = outlet_new(&x->x_ob, gensym("list"));
    return (void *)x;
}

//float foo(float i) {return 20*pow(2,i/3);}
//float unf(float f) {return log(f/20)/log(2)*3;}
static void iso226_float(t_iso226 *x, t_float ph) {
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

    if ((ph < 0) || (ph > 90)) {pd_error(x,"Phon value out of bounds!"); return;}
    //Deriving sound pressure level from loudness level (iso226 sect 4.1)
    float Af[29];
    for (i=0; i<29; i++) Af[i] = 4.47E-3 * (pow(10,0.025*ph) - 1.15) + pow(0.4*pow(10,((Tf[i]+Lu[i])/10)-9 ),af[i]);
    for (i=0; i<29; i++) SETFLOAT(&spl[i], ((10/af[i]) * log10(Af[i])) - Lu[i] + 94);

    outlet_list(x->x_outlet,&s_list,29,spl);
}
static double phons[11]={2,10,20,30,40,50,60,70,80,90,100};
static int eqlbandbins[43]= {
      1,  2,  3,  4,  5,  6,  7,  8,  9, 11, 13, 15, 17, 19,
     22, 25, 28, 32, 36, 41, 46, 52, 58, 65, 73, 82, 92,103,
    116,129,144,161,180,201,225,251,280,312,348,388,433,483,513};
static float logfreqs[43];
static float contours[42][11]= {
{ 47.88, 59.68, 68.55, 75.48, 81.71, 87.54, 93.24, 98.84,104.44,109.94,115.31}, // 0
{ 29.04, 41.78, 51.98, 60.18, 67.51, 74.54, 81.34, 87.97, 94.61,101.21,107.74},
{ 20.72, 32.83, 43.44, 52.18, 60.24, 67.89, 75.34, 82.70, 89.97, 97.23,104.49},
{ 15.87, 27.14, 37.84, 46.94, 55.44, 63.57, 71.51, 79.34, 87.14, 94.97,102.37},
{ 12.64, 23.24, 33.91, 43.27, 52.07, 60.57, 68.87, 77.10, 85.24, 93.44,100.90},
{ 10.31, 20.43, 31.03, 40.54, 49.59, 58.33, 66.89, 75.43, 83.89, 92.34,100.80}, // 5
{  8.51, 18.23, 28.83, 38.41, 47.65, 56.59, 65.42, 74.16, 82.89, 91.61,100.33},
{  7.14, 16.55, 27.11, 36.79, 46.16, 55.27, 64.29, 73.24, 82.15, 91.06, 99.97},
{  5.52, 14.58, 25.07, 34.88, 44.40, 53.73, 62.95, 72.18, 81.31, 90.44, 99.57},
{  3.98, 12.69, 23.10, 32.99, 42.69, 52.27, 61.66, 71.15, 80.54, 89.93, 99.31},
{  2.99, 11.43, 21.76, 31.73, 41.49, 51.22, 60.88, 70.51, 80.11, 89.70, 99.30}, // 10
{  2.35, 10.58, 20.83, 30.86, 40.68, 50.51, 60.33, 70.08, 79.83, 89.58, 99.32},
{  2.05, 10.12, 20.27, 30.35, 40.22, 50.10, 59.97, 69.82, 79.67, 89.52, 99.38},
{  2.00,  9.93, 20.00, 30.07, 40.00, 49.93, 59.87, 69.80, 79.73, 89.67, 99.60},
{  2.19, 10.00, 20.00, 30.00, 40.00, 50.00, 59.99, 69.99, 79.98, 89.98, 99.97}, // 14 : 1 KHz
{  2.71, 10.56, 20.61, 30.71, 40.76, 50.81, 60.86, 70.96, 81.01, 91.06,101.17},
{  3.11, 11.05, 21.19, 31.41, 41.53, 51.64, 61.75, 71.95, 82.05, 92.15,102.33},
{  2.39, 10.69, 21.14, 31.52, 41.73, 51.95, 62.11, 72.31, 82.46, 92.56,102.59},
{  1.50, 10.11, 20.82, 31.32, 41.62, 51.92, 62.12, 72.32, 82.52, 92.63,102.56},
{ -0.17,  8.50, 19.27, 29.77, 40.07, 50.37, 60.57, 70.77, 80.97, 91.13,101.23},
{ -1.80,  6.96, 17.77, 28.29, 38.61, 48.91, 59.13, 69.33, 79.53, 89.71, 99.86},
{ -3.42,  5.49, 16.36, 26.94, 37.31, 47.61, 57.88, 68.08, 78.28, 88.41, 98.39},
{ -4.73,  4.38, 15.34, 25.99, 36.39, 46.71, 57.01, 67.21, 77.41, 87.51, 97.41},
{ -5.73,  3.63, 14.74, 25.48, 35.88, 46.26, 56.56, 66.76, 76.96, 87.06, 96.96},
{ -6.24,  3.33, 14.59, 25.39, 35.84, 46.22, 56.52, 66.72, 76.92, 87.04, 97.00},
{ -6.09,  3.62, 15.03, 25.83, 36.37, 46.70, 57.00, 67.20, 77.40, 87.57, 97.68},
{ -5.32,  4.44, 15.90, 26.70, 37.28, 47.60, 57.90, 68.10, 78.30, 88.52, 98.78},
{ -3.49,  6.17, 17.52, 28.32, 38.85, 49.22, 59.52, 69.72, 79.92, 90.20,100.61},
{ -0.81,  8.58, 19.73, 30.44, 40.90, 51.24, 61.52, 71.69, 81.87, 92.15,102.63},
{  2.91, 11.82, 22.64, 33.17, 43.53, 53.73, 63.96, 74.09, 84.22, 94.45,104.89},
{  6.68, 15.19, 25.71, 36.03, 46.25, 56.31, 66.45, 76.49, 86.54, 96.72,107.15},
{ 10.43, 18.65, 28.94, 39.02, 49.01, 58.98, 68.93, 78.78, 88.69, 98.83,109.36},
{ 13.56, 21.65, 31.78, 41.68, 51.45, 61.31, 71.07, 80.73, 90.48,100.51,111.01},
{ 14.36, 22.91, 33.19, 43.09, 52.71, 62.37, 71.92, 81.38, 90.88,100.56,110.56},
{ 15.06, 23.90, 34.23, 44.05, 53.48, 62.90, 72.21, 81.43, 90.65, 99.93,109.34},
{ 15.36, 23.90, 33.89, 43.31, 52.40, 61.42, 70.29, 79.18, 88.00, 96.69,105.17},
{ 15.60, 23.90, 33.60, 42.70, 51.50, 60.20, 68.70, 77.30, 85.80, 94.00,101.70},
{ 15.60, 23.90, 33.60, 42.70, 51.50, 60.20, 68.70, 77.30, 85.80, 94.00,101.70},
{ 15.60, 23.90, 33.60, 42.70, 51.50, 60.20, 68.70, 77.30, 85.80, 94.00,101.70},
{ 15.60, 23.90, 33.60, 42.70, 51.50, 60.20, 68.70, 77.30, 85.80, 94.00,101.70},
{ 15.60, 23.90, 33.60, 42.70, 51.50, 60.20, 68.70, 77.30, 85.80, 94.00,101.70},
{ 15.60, 23.90, 33.60, 42.70, 51.50, 60.20, 68.70, 77.30, 85.80, 94.00,101.70}};

static void iso226b_left(t_func *x, t_float left) {
        float freq = left;
        if (freq > 12500) {freq = 12500;}
        if (freq < 20) {freq = 20;}
        freq = (log(freq) - log(20)) / (log(12500) - log(20)); // make freq a number between 0 and 1
		float lastfreq = logfreqs[0];
        int i;
        float fprop=0;
        for (i=1; i<43; i++) {
            float val = logfreqs[i];
            if(freq >= lastfreq && freq <= val) {
                fprop = (freq - lastfreq) / (val - lastfreq);
                break;
            }
			lastfreq = val;
		}
		float contour[11];
		int j;
		for (j=0; j<11; j++) contour[j] = (1-fprop) * contours[i-1][j] + fprop * contours[i][j];

        float db = x->right;
	    if(db<contour[0]) db=0;
		else if (db>contour[10]) db=phons[10];
		else {
			float prop=0.0;
			for (j=1; j<11; ++j) {
				if(db<contour[j]) {
					prop= (db-contour[j-1])/(contour[j]-contour[j-1]);
					break;
				}
				if(j==10) prop=1.0;
			}
			 outlet_float(x->x_outlet,(1.f-prop)*phons[j-1]+ prop*phons[j]);
			//printf("prop %f db %f j %d\n",prop,db,j);
		}
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%

void psycho_setup(void){
    roughness_class = class_new(gensym("roughness"), (t_newmethod)roughness_new, (t_method)roughness_free, sizeof(t_roughness), 0, A_GIMME, 0);
	class_addlist(roughness_class, roughness_freqs);
	class_addmethod(roughness_class, (t_method)roughness_amps, gensym("amps"), A_GIMME, 0); 
    class_addmethod(roughness_class, (t_method)roughness_phonsteps, gensym("phonsteps"), A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_loudness,  gensym("loudness"),  A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_mean,      gensym("weight"),    A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_curve,     gensym("curve"),     A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_masking,   gensym("masking"),   A_FLOAT, 0);
    class_addmethod(roughness_class, (t_method)roughness_masking,   gensym("sones"),     A_FLOAT, 0);
	
	flunson_class = class_new(gensym("flunson"), (t_newmethod)flunson_new,
                              (t_method)flunson_free, sizeof(t_flunson), 0, 0);
	class_addlist(flunson_class, flunson_freqs);
	class_addanything(flunson_class, flunson_amps); // use class_addmethod with gensym("amps")
    
    phon_class = class_new(gensym("phon"), (t_newmethod)phon_new, 0, sizeof(t_func), 0, A_DEFFLOAT, 0);
    class_addmethod(phon_class,(t_method)function_right,gensym("right"),A_FLOAT,0);
    class_addfloat(phon_class,(t_method)phon_left);
    
    dbA_class = class_new(gensym("dbA"),(t_newmethod)dbA_new,0,sizeof(t_func),0,A_DEFFLOAT,0);
    class_addmethod(dbA_class,(t_method)function_right,gensym("right"),A_FLOAT,0);
    class_addfloat(dbA_class,(t_method)dbA_left);
    
    iso226b_class = class_new(gensym("iso226b"),(t_newmethod)iso226b_new,0,sizeof(t_func),0,A_DEFFLOAT,0);
    class_addmethod(iso226b_class,(t_method)function_right,gensym("right"),A_FLOAT,0);
    class_addfloat(iso226b_class,(t_method)iso226b_left);
    
    iso226_class = class_new(gensym("iso226"),(t_newmethod)iso226_new,0,sizeof(t_func),0, 0);
    class_addfloat(iso226_class,(t_method)iso226_float);
    
    post("Psycho library v0.0.1 alpha-0 (unreleased)");
    for(int i = 0; i < 43; i++)
        logfreqs[i] = log(eqlbandbins[i]/512.0*44100.0); // ????????????
}
