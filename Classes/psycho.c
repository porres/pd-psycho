
#include <string.h>
#include <math.h>
#include "m_pd.h"

//!!!!! PSYCHOACOUSTICAL FUNCTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!/

//=============================== Amplitude & Loudness ====================================/

//***** Convert Phons/hz to dBA (Barlow's function, Robson/Dadson's data)
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

// !!!!!!!!!!!!!!!!!!!!!! OBJECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//************************ FLUNSON OBJECT ************************/

t_class *flunson_class;

typedef struct flunson{
    t_object x_ob;
    t_outlet *x_outlet;
	int ac;
	float *amps;
}t_flunson;

// ******* [rmstodb] - Converts Linear Amplitude to dB(A)
static float gain2db(float amp){
    if(amp > 0.00001)
        return 100+20*log10f(amp);
    return 0;
}

// ******* [dbtorms] - Converts dB(A) to Linear Amplitude
static float db2gain(float db){
    return pow(10, (db-100)/20);
}

void flunson_freqs(t_flunson *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(ac > x->ac)
        ac = x->ac; // too many args.
    t_atom outs[ac];
    int i;
//  for (i=0; i<ac; i++) SETFLOAT(outs+i, phon(atom_getfloat(av+i), x->amps[i] ) );
    for(i = 0; i < ac; i++)
        SETFLOAT(outs+i, db2gain(phon(atom_getfloat(av+i), gain2db(x->amps[i]), 20)));
//  for (i=0; i<ac; i++) SETFLOAT(outs+i,db2gain( dbA(atom_getfloat(av+i),gain2db(x->amps[i]))));
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

//************************ Phon / dBA / iso226b ************************/
typedef struct func{
    t_object x_ob;
    t_outlet *x_outlet;
    float right;
}t_func;

static t_class *db2phon_class;

static void *db2phon_new(t_float right){
    t_func *x = (t_func *)pd_new(db2phon_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("float"), gensym("right"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->right = right;
    return (void *)x;
}

static t_class *phon2db_class;

static void *phon2db_new(t_float right){
    t_func *x = (t_func *)pd_new(phon2db_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("float"), gensym("right"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->right = right;
    return (void *)x;
}

static t_class *iso226b_class;

static void *iso226b_new(t_float right){
    t_func *x = (t_func *)pd_new(iso226b_class);
    inlet_new(&x->x_ob, &x->x_ob.ob_pd, gensym("float"), gensym("right"));
    x->x_outlet = outlet_new(&x->x_ob, gensym("float"));
    x->right = right;
    return (void *)x;
}

static void function_right (t_func *x, float right){
    x->right = right;
}

static void phon_left (t_func *x, float left){
    outlet_float(x->x_outlet, phon(left, x->right, 20));
}

static void dbA_left (t_func *x, float left){
    outlet_float(x->x_outlet, dbA(left, x->right));
}

static double phons[11] = {2,10,20,30,40,50,60,70,80,90,100};

static int eqlbandbins[43]= {
      1,  2,  3,  4,  5,  6,  7,  8,  9, 11, 13, 15, 17, 19,
     22, 25, 28, 32, 36, 41, 46, 52, 58, 65, 73, 82, 92,103,
    116,129,144,161,180,201,225,251,280,312,348,388,433,483,513
};

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

static void iso226b_left(t_func *x, t_float left){
        float freq = left;
        if (freq > 12500) {freq = 12500;}
        if (freq < 20) {freq = 20;}
        freq = (log(freq) - log(20)) / (log(12500) - log(20)); // make freq a number between 0 and 1
		float lastfreq = logfreqs[0];
        int i;
        float fprop=0;
        for(i = 1; i < 43; i++){
            float val = logfreqs[i];
            if(freq >= lastfreq && freq <= val){
                fprop = (freq - lastfreq) / (val - lastfreq);
                break;
            }
			lastfreq = val;
		}
		float contour[11];
		int j;
		for(j = 0; j < 11; j++) contour[j] = (1-fprop) * contours[i-1][j] + fprop * contours[i][j];
        float db = x->right;
	    if(db < contour[0])
            db = 0;
		else if (db > contour[10])
            db = phons[10];
		else{
			float prop = 0.0;
			for(j = 1; j < 11; ++j){
				if(db < contour[j]){
					prop= (db-contour[j-1])/(contour[j]-contour[j-1]);
					break;
				}
				if(j == 10)
                    prop = 1.0;
			}
            outlet_float(x->x_outlet, (1.f-prop) * phons[j-1] + prop*phons[j]);
		}
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%
void psycho_setup(void){
	flunson_class = class_new(gensym("flunson"), (t_newmethod)flunson_new,
                (t_method)flunson_free, sizeof(t_flunson), 0, 0);
	class_addlist(flunson_class, flunson_freqs);
	class_addanything(flunson_class, flunson_amps); // use class_addmethod with gensym("amps")
    
    db2phon_class = class_new(gensym("db2phon"), (t_newmethod)db2phon_new, 0, sizeof(t_func), 0, A_DEFFLOAT, 0);
    class_addmethod(db2phon_class,(t_method)function_right,gensym("right"), A_FLOAT, 0);
    class_addfloat(db2phon_class,(t_method)phon_left);
    
    phon2db_class = class_new(gensym("phon2db"),(t_newmethod)phon2db_new,0,sizeof(t_func), 0, A_DEFFLOAT, 0);
    class_addmethod(phon2db_class, (t_method)function_right, gensym("right"), A_FLOAT, 0);
    class_addfloat(phon2db_class, (t_method)dbA_left);
    
    iso226b_class = class_new(gensym("iso226b"), (t_newmethod)iso226b_new, 0,
                              sizeof(t_func), 0, A_DEFFLOAT, 0);
    class_addmethod(iso226b_class, (t_method)function_right, gensym("right"), A_FLOAT, 0);
    class_addfloat(iso226b_class,(t_method)iso226b_left);
    
    for(int i = 0; i < 43; i++)
        logfreqs[i] = log(eqlbandbins[i]/512.0*44100.0); // ????????????
    
    post("Psycho library v0.0.1 alpha-0 (unreleased)");
}
