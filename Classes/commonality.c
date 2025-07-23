/*
 By porres, 2010-2019, based on:
 Parncutt, Richard & Strasburger, Hans. (1994).
 Applying psychoacoustics in composition: Harmonic progressions of non-harmonic sonorities.
 Perspectives of New Music. 32. 88-129. 10.2307/833600.
*/

#include "m_pd.h"
#include <math.h>

typedef struct floatlist{
	int         n;
	t_float    *v;
}t_floatlist;

//// proxy
typedef struct proxy_obj{
	t_pd        x_obj;
	t_floatlist list;
}t_proxy_obj;

t_class *proxy_obj_class;

void proxy_obj_list(t_proxy_obj *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
    if(x->list.v)
        freebytes(x->list.v, x->list.n*sizeof(float));
    x->list.v = (float *)getbytes(ac*sizeof(float));
    x->list.n = ac;
    for(int i = 0; i < ac; i++)
        x->list.v[i] = atom_getfloat(av+i);
}

void proxy_obj_free(t_proxy_obj *x){
    freebytes(x->list.v, x->list.n*sizeof(float));
}

void *proxy_obj_new(void){
    t_proxy_obj *x = (t_proxy_obj *)pd_new(proxy_obj_class);
    x->list.n = 0;
    x->list.v = 0;
    return(x);
}

// Pitch Commonality ---------------------------------------------------------------------------

typedef struct commonality{
	t_object        x_obj;
	t_outlet       *outlet0;
	t_proxy_obj   **proxies;    // one spectrum per proxy
	int             n;          // number of spectra (and of proxies)
}t_commonality;

t_class *commonality_class;

// x and y are arrays with m elements / corr is a single variable.
void pearson(int m, float *x, float *y, float *corr){
    float sumx = 0, sumy = 0, sumxx = 0, sumyy = 0, sumxy = 0, sumrootxy = 0;
    // omits bottom and top half octaves of range of hearing, otherwise as i=0 gives infinite values (bug!)
//    for(int i = 6; i <= m-6; i++){
    for(int i = 0; i < m; i++){
        sumx += x[i];
        if(x[i] < 0)
            post("[commonality]: negative value: i=%d x=%5.0f", i, x[i]);
        sumy += y[i];
        if(y[i] < 0)
            post("[commonality]: negative value: i=%d y=%5.0f", i, y[i]);
        sumxx += x[i]*x[i];
        sumyy += y[i]*y[i];
        sumxy += x[i]*y[i];
        sumrootxy += pow(x[i]*y[i], 0.5);
    }
    float sx = sqrt((sumxx - sumx*sumx / m) / (m-1)); // stdev(x)
    float sy = sqrt((sumyy - sumy*sumy / m) / (m-1)); // stdev(y)
    float sxy =     (sumxy - sumx*sumy / m) / (m-1);  // covariance(x, y)
    *corr = sxy / (sx*sy); // pearson correlation coefficient r */
}

void commonality_bang(t_commonality *x){
	int n = x->n, m = x->proxies[0]->list.n;
	int i, j;
	for(i = 1; i < n; i++)
        if(m != x->proxies[i]->list.n){
        pd_error(x, "[commonality]: input lists need to have the same");
		return;
	}
	float pitCorr[10][10]; // pitch correlation
	for(i = 0; i < 10; i++)
        pitCorr[i][i] = 1;
	// PART 3. CALCULATE PITCH RELATIONSHIPS BETWEEN SONORITIES.
	for(i = 0; i < n; i++) { // symmetric matrix, so we do it for a triangular matrix and copy to the other half
		for(j = i+1; j < n; j++){
			float corr;
			pearson(m, x->proxies[i]->list.v,
                        x->proxies[j]->list.v, &corr);
			pitCorr[i][j] = pitCorr[j][i] = corr;
		}
	}
	if(n == 2){
		outlet_float(x->outlet0, pitCorr[0][1]);
	}
    else{
		t_atom outlist[n*n];
		for(i = 0; i < n; i++)
            for(j = 0; j < n; j++)
                SETFLOAT(outlist+i*n+j, pitCorr[i][j]);
		outlet_list(x->outlet0, &s_list, n*n, outlist);
	}
}

void commonality_list(t_commonality *x, t_symbol *s, int ac, t_atom *av){
    t_symbol *dummy = s;
    dummy = NULL;
	proxy_obj_list(x->proxies[0], &s_list, ac, av);
	commonality_bang(x);
}

void commonality_free(t_commonality *x){
    for(int i = 0; i < x->n; i++)
        pd_free((t_pd *)x->proxies[i]);
}

void *commonality_new(t_floatarg f){
    t_commonality *x = (t_commonality *)pd_new(commonality_class);
    x->n = f < 2 ? 2 : f > 10 ? 10 : f;
    x->proxies = (t_proxy_obj **)getbytes(x->n*sizeof(t_proxy_obj *));
    for(int i = 0; i < x->n; i++){
        x->proxies[i] = proxy_obj_new();
        if(i)
            inlet_new((t_object *)x, (t_pd *)x->proxies[i], &s_list, &s_list);
    }
    x->outlet0 = outlet_new((t_object *)x, &s_float);
    return (void *)x;
}

void commonality_setup(void){
    commonality_class = class_new(gensym("commonality"), (t_newmethod)commonality_new,
            (t_method)commonality_free, sizeof(t_commonality), 0, A_DEFFLOAT, 0);
    class_addlist(commonality_class, (t_method)commonality_list);
    proxy_obj_class = class_new(gensym("proxy_obj"), (t_newmethod)proxy_obj_new,
            (t_method)proxy_obj_free, sizeof(t_proxy_obj), CLASS_PD, 0);
	class_addlist(proxy_obj_class, (t_method)proxy_obj_list);
}
