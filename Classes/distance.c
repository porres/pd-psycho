/*
 By porres, 2010-2019, based on:
 Parncutt, Richard & Strasburger, Hans. (1994).
 Applying psychoacoustics in composition: Harmonic progressions of non-harmonic sonorities.
 Perspectives of New Music. 32. 88-129. 10.2307/833600.
 */

#include "m_pd.h"
#include <math.h>
#include <stdlib.h>

typedef struct floatlist{
	int n;
	t_float *v;
}floatlist;

// proxy

typedef struct proxy_obj{
    t_pd o;
    floatlist list;
}proxy_obj;

t_class *proxy_obj_class;

void proxy_obj_list(proxy_obj *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->list.v)
        freebytes(x->list.v, x->list.n*sizeof(float));
    x->list.v = (float *)getbytes(ac*sizeof(float));
    x->list.n = ac;
    for(int i = 0; i < ac; i++)
        x->list.v[i] = atom_getfloat(av+i);
}

void proxy_obj_free (proxy_obj *x){
    freebytes(x->list.v, x->list.n*sizeof(float));
}

void *proxy_obj_new(void){
    proxy_obj *x = (proxy_obj *)pd_new(proxy_obj_class);
    x->list.n = 0;
    x->list.v = 0;
    return(x);
}

// Pitch Distance ------------------------------------------------------------------------------

typedef struct distance{
	t_object o;
	t_outlet *outlet0;
	proxy_obj **proxies; // one spectrum per proxy
	int n; // number of spectra (and of proxies)
}distance;

t_class *distance_class;

// a and b are two arrays of size m
float pitch_distance (int m, float *a, float *b){
	float sumab = 0,sumaa = 0, sumbb = 0;
	for(int i = 0; i < m; i++){
		for(int j = 0; j < m; j++){
			sumab += a[i]*b[j]*abs(i-j);
			sumaa += a[i]*a[j]*abs(i-j);
			sumbb += b[i]*b[j]*abs(i-j);
		}
	}
	return (sumab - sqrt(sumaa*sumbb)) * 120/m;
}

// D = sum(i,sum(j,a[i]b[j]|j-i])) - sqrt(sum(i,sum(j,a[i]a[j]|j-i|)) sum(i,sum(j,b[i]b[j]|j-i|))
void distance_bang (distance *x) {
	int n = x->n, m = x->proxies[0]->list.n;
	int i, j;
	for(i = 1; i < n; i++)
        if(m != x->proxies[i]->list.n){
		pd_error(x,"list of inlet %d does not have the same size (%d) as list of inlet 0 (%d)",
            i, m, x->proxies[i]->list.n);
		return;
	}
	float pitDist[10][10]; // pitch distance
	// (d) initialize pitCom array
	for(i = 0; i < 10; i++)
        pitDist[i][i] = 0;
	// PART 3. CALCULATE PITCH RELATIONSHIPS BETWEEN SONORITIES.
	for(i = 0; i < n; i++){ // symmetric matrix, so we do it for a triangular matrix and copy to the other half
		for(j = i+1; j < n; j++){
			pitDist[i][j] = pitDist[j][i] = pitch_distance(m,
                x->proxies[i]->list.v,x->proxies[j]->list.v);
		}
	}
	if(n == 2)
		outlet_float(x->outlet0,pitDist[0][1]);
    else{
		t_atom outlist[n*n];
		for(i = 0; i < n; i++)
            for(j = 0; j < n; j++)
                SETFLOAT(outlist+i*n+j, pitDist[i][j]);
		outlet_list(x->outlet0, &s_list, n*n, outlist);
	}
}

void distance_list(distance *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
	proxy_obj_list(x->proxies[0], &s_list, ac, av);
	distance_bang(x);
}

void distance_free(distance *x){
    for(int i = 0; i < x->n; i++)
        pd_free((t_pd *)x->proxies[i]);
}

void *distance_new(t_floatarg f){
    distance *x = (distance *)pd_new(distance_class);
    int n = f < 2 ? 2 : f > 10 ? 10 : f;
    x->n = n;
    x->proxies = (proxy_obj **)getbytes(n*sizeof(proxy_obj *));
    for(int i = 0; i < n; i++){
        x->proxies[i] = proxy_obj_new();
        if(i)
            inlet_new((t_object *)x, (t_pd *)x->proxies[i], &s_list, &s_list);
    }
    x->outlet0 = outlet_new((t_object *)x,gensym("list"));
    return(x);
}

// setup

void distance_setup(void){ // (MAXTA[cat]) -> (SalTA, mul)
    proxy_obj_class = class_new(gensym("proxy_obj"), (t_newmethod)proxy_obj_new,
            (t_method)proxy_obj_free, sizeof(proxy_obj), CLASS_PD, 0);
	class_addlist(proxy_obj_class, (t_method)proxy_obj_list);
    distance_class = class_new(gensym("distance"), (t_newmethod)distance_new,
            (t_method)distance_free, sizeof(distance), 0, A_DEFFLOAT, 0);
	class_addlist(distance_class, (t_method)distance_list);
}
