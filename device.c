/* c-ray-mt - a simple multithreaded raytracing filter.
 * Copyright (C) 2006 John Tsiombikas <nuclear@siggraph.org>
 *
 * You are free to use, modify and redistribute this program under the
 * terms of the GNU General Public License v2 or (at your option) later.
 * see "http://www.gnu.org/licenses/gpl.txt" for details.
 */

/*
 * Modified to support the Adapteva Epiphany architecure using the COPRTHR-2
 * interface.
 * Copyright (C) 2016 James A. Ross (james.a.ross@gmail.com)
 */

#include <coprthr2.h>
#include "common.h"

#define RAY_MAG        1000.0f       // trace rays of this magnitude
#define MAX_RAY_DEPTH  5             // raytrace recursion limit
#define FOV            0.78539816f   // field of view in rads (pi/4)
#define INV_HALF_FOV   2.546479089f
#define ERR_MARGIN     1.0e-5f       // an arbitrary error margin to avoid surface acne
#define M_LN2          0.6931471805f
#define M_INVLN2       1.442695041f  // 1/M_LN2

/* some helpful macros... */
#define SQ(x) ((x) * (x))
#define MAX(a, b)	((a) > (b) ? (a) : (b))
#define MIN(a, b)	((a) < (b) ? (a) : (b))
#define DOT(a, b)	((a).x * (b).x + (a).y * (b).y + (a).z * (b).z)
#define NORMALIZE(a)	do {\
	float ilen = _rsqrtf(DOT(a, a)); \
	(a).x *= ilen; (a).y *= ilen; (a).z *= ilen;\
} while(0);

void render_scanline(int xsz, int ysz, int sl, uint32_t *fb, int samples);
static struct vec3 trace(struct ray ray, int depth);
static struct vec3 shade(struct sphere *obj, struct spoint *sp, int depth);
static struct vec3 reflect(struct vec3 v, struct vec3 n);
static struct vec3 cross_product(struct vec3 v1, struct vec3 v2);
static struct ray get_primary_ray(int x, int y, int sample);
static struct vec3 get_sample_pos(int x, int y, int sample);
static struct vec3 jitter(int x, int y, int s);
static int ray_sphere(const struct sphere *sph, struct ray ray, struct spoint *sp);

/* Similar newton iteration schemes are used for sqrt, rsqrt, and inv routines */
static inline float _sqrtf( const float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5f;
	x2 = number * 0.5f;
	y	= number;
	i	= * ( long * ) &y;
	i	= 0x5f3759df - ( i >> 1 );
	y	= * ( float * ) &i;
	y	= y * ( threehalfs - ( x2 * y * y ) );
	y	= y * ( threehalfs - ( x2 * y * y ) );
	y	= y * ( threehalfs - ( x2 * y * y ) );
	y	= y * ( threehalfs - ( x2 * y * y ) );
	return y * number;
}

static inline float _rsqrtf( const float number )
{
	long i;
	float x2, y;
	const float threehalfs = 1.5f;
	x2 = number * 0.5f;
	y	= number;
	i	= * ( long * ) &y;
	i	= 0x5f3759df - ( i >> 1 );
	y	= * ( float * ) &i;
	y	= y * ( threehalfs - ( x2 * y * y ) );
	y	= y * ( threehalfs - ( x2 * y * y ) );
	y	= y * ( threehalfs - ( x2 * y * y ) );
	y	= y * ( threehalfs - ( x2 * y * y ) );
	return y;
}

static inline float __inv(const float number)
{
	union {
		float f;
		uint32_t x;
	} u = {number};
	u.x = 0x7EEEEBB3 - u.x;
	u.f = u.f * (2.0f - u.f * number);
	u.f = u.f * (2.0f - u.f * number);
	u.f = u.f * (2.0f - u.f * number);
	u.f = u.f * (2.0f - u.f * number);
	return u.f;
}

/*
 * Copyright applies for copysignf, scalbnf, powf, and supporting code:
 */
/*
 * ===========================================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.  Permission to use,
 * copy, modify, and distribute this software is freely granted, provided that
 * this notice is preserved.
 * ===========================================================================
 */

typedef union
{
	float value;
	uint32_t word;
} ieee_float_shape_type;

#define GET_FLOAT_WORD(i,d) \
do { \
	ieee_float_shape_type gf_u; \
	gf_u.value = (d); \
	(i) = gf_u.word; \
} while (0)

#define SET_FLOAT_WORD(d,i) \
do { \
	ieee_float_shape_type sf_u; \
	sf_u.word = (i); \
	(d) = sf_u.value; \
} while (0)

static const float
bp[]    = {1.0, 1.5,},
dp_h[]  = { 0.0, 5.84960938e-01,}, // 0x3f15c000
dp_l[]  = { 0.0, 1.56322085e-06,}, // 0x35d1cfdc
zero	  = 0.0,
one	  = 1.0,
two	  = 2.0,
two24	  = 16777216.0,        // 0x4b800000
huge	  = 1.0e30,
tiny	  = 1.0e-30,
two25   = 3.355443200e+07,	  // 0x4c000000
twom25  = 2.9802322388e-08,  // 0x33000000
// poly coefs for (3/2)*(log(x)-2s-2/3*s**3
L1      =  6.0000002384e-01, // 0x3f19999a
L2      =  4.2857143283e-01, // 0x3edb6db7
L3      =  3.3333334327e-01, // 0x3eaaaaab
L4      =  2.7272811532e-01, // 0x3e8ba305
L5      =  2.3066075146e-01, // 0x3e6c3255
L6      =  2.0697501302e-01, // 0x3e53f142
P1      =  1.6666667163e-01, // 0x3e2aaaab
P2      = -2.7777778450e-03, // 0xbb360b61
P3      =  6.6137559770e-05, // 0x388ab355
P4      = -1.6533901999e-06, // 0xb5ddea0e
P5      =  4.1381369442e-08, // 0x3331bb4c
lg2     =  6.9314718246e-01, // 0x3f317218
lg2_h   =  6.93145752e-01,   // 0x3f317200
lg2_l   =  1.42860654e-06,   // 0x35bfbe8c
ovt     =  4.2995665694e-08, // -(128-log2(ovfl+.5ulp))
cp	     =  9.6179670095e-01, // 0x3f76384f =2/(3ln2)
cp_h    =  9.6191406250e-01, // 0x3f764000 =12b cp
cp_l    = -1.1736857402e-04, // 0xb8f623c6 =tail of cp_h
ivln2	  =  1.4426950216e+00, // 0x3fb8aa3b =1/ln2
ivln2_h =  1.4426879883e+00, // 0x3fb8aa00 =16b 1/ln2
ivln2_l =  7.0526075433e-06; // 0x36eca570 =1/ln2 tail

static float
__copysignf(float x, float y)
{
	uint32_t ix,iy;
	GET_FLOAT_WORD(ix,x);
	GET_FLOAT_WORD(iy,y);
	SET_FLOAT_WORD(x,(ix&0x7fffffff)|(iy&0x80000000));
	return x;
}

static inline float
__scalbnf(float x, int n)
{
	int32_t k,ix;
	GET_FLOAT_WORD(ix,x);
	k = (ix&0x7f800000)>>23; /* extract exponent */
	if (k==0) { /* 0 or subnormal x */
		if ((ix&0x7fffffff)==0) return x; /* +-0 */
		x *= two25;
		GET_FLOAT_WORD(ix,x);
		k = ((ix&0x7f800000)>>23) - 25;
		if (n< -50000) return tiny*x; /*underflow*/
	}
	if (k==0xff) return x+x; /* NaN or Inf */
	k = k+n;
	if (k > 0xfe) return huge*__copysignf(huge,x); /* overflow */
	if (k > 0) /* normal result */
	{SET_FLOAT_WORD(x,(ix&0x807fffff)|(k<<23)); return x;}
	if (k <= -25) {
		if (n > 50000) /* in case integer overflow in n+k */
			return huge*__copysignf(huge,x); /*overflow*/
		else return tiny*__copysignf(tiny,x); /*underflow*/
	}
	k += 25; /* subnormal result */
	SET_FLOAT_WORD(x,(ix&0x807fffff)|(k<<23));
	return x*twom25;
}

static inline float
__powf(float x, float y)
{
	float z,ax,z_h,z_l,p_h,p_l;
	float y1,t1,t2,r,s,sn,t,u,v,w;
	int32_t i,j,k,yisint,n;
	int32_t hx,hy,ix,iy,is;

	GET_FLOAT_WORD(hx,x);
	GET_FLOAT_WORD(hy,y);
	ix = hx&0x7fffffff;
	iy = hy&0x7fffffff;

	/* y==zero: x**0 = 1 */
	if(iy==0) return one;

	/* x==1: 1**y = 1, even if y is NaN */
	if (hx==0x3f800000) return one;

	/* y!=zero: result is NaN if either arg is NaN */
	if(ix > 0x7f800000 || iy > 0x7f800000)
		return (x+0.0F)+(y+0.0F);

	/* determine if y is an odd int when x < 0
	 * yisint = 0	... y is not an integer
	 * yisint = 1	... y is an odd int
	 * yisint = 2	... y is an even int
	 */
	yisint = 0;
	if(hx<0) {
		if(iy>=0x4b800000) yisint = 2; /* even integer y */
		else if(iy>=0x3f800000) {
		k = (iy>>23)-0x7f; /* exponent */
		j = iy>>(23-k);
		if((j<<(23-k))==iy) yisint = 2-(j&1);
		}
	}

	ax = *(float*)&ix;

	n = ((uint32_t)hx>>31)-1;

	/* (x<0)**(non-int) is NaN */
	if((n|yisint)==0) return 0.0f/0.0f;

	sn = one; /* s (sign of result -ve**odd) = -1 else = 1 */
	if((n|(yisint-1))==0) sn = -one;/* (-ve)**(odd int) */
	{
		float s2,s_h,s_l,t_h,t_l;
		n = 0;
		/* take care subnormal number */
		if(ix<0x00800000)
		{ax *= two24; n -= 24; GET_FLOAT_WORD(ix,ax); }
		n += ((ix)>>23)-0x7f;
		j = ix&0x007fffff;
		/* determine interval */
		ix = j|0x3f800000;		/* normalize ix */
		if(j<=0x1cc471) k=0;	/* |x|<sqrt(3/2) */
		else if(j<0x5db3d7) k=1;	/* |x|<sqrt(3) */
		else {k=0;n+=1;ix -= 0x00800000;}
		SET_FLOAT_WORD(ax,ix);

		/* compute s = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
		u = ax-bp[k];		/* bp[0]=1.0, bp[1]=1.5 */
		v = __inv(ax+bp[k]);
		s = u*v;
		s_h = s;
		GET_FLOAT_WORD(is,s_h);
		SET_FLOAT_WORD(s_h,is&0xfffff000);
		/* t_h=ax+bp[k] High */
		is = ((ix>>1)&0xfffff000)|0x20000000;
		SET_FLOAT_WORD(t_h,is+0x00400000+(k<<21));
		t_l = ax - (t_h-bp[k]);
		s_l = v*((u-s_h*t_h)-s_h*t_l);
		/* compute log(ax) */
		s2 = s*s;
		r = s2*s2*(L1+s2*(L2+s2*(L3+s2*(L4+s2*(L5+s2*L6)))));
		r += s_l*(s_h+s);
		s2 = s_h*s_h;
		t_h = 3.0f+s2+r;
		GET_FLOAT_WORD(is,t_h);
		SET_FLOAT_WORD(t_h,is&0xfffff000);
		t_l = r-((t_h-3.0f)-s2);
		/* u+v = s*(1+...) */
		u = s_h*t_h;
		v = s_l*t_h+t_l*s;
		/* 2/(3log2)*(s+...) */
		p_h = u+v;
		GET_FLOAT_WORD(is,p_h);
		SET_FLOAT_WORD(p_h,is&0xfffff000);
		p_l = v-(p_h-u);
		z_h = cp_h*p_h;		/* cp_h+cp_l = 2/(3*log2) */
		z_l = cp_l*p_h+p_l*cp+dp_l[k];
		/* log2(ax) = (s+..)*2/(3*log2) = n + dp_h + z_h + z_l */
		t = (float)n;
		t1 = (((z_h+z_l)+dp_h[k])+t);
		GET_FLOAT_WORD(is,t1);
		SET_FLOAT_WORD(t1,is&0xfffff000);
		t2 = z_l-(((t1-t)-dp_h[k])-z_h);
	}

	/* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
	GET_FLOAT_WORD(is,y);
	SET_FLOAT_WORD(y1,is&0xfffff000);
	p_l = (y-y1)*t1+y*t2;
	p_h = y1*t1;
	z = p_l+p_h;
	GET_FLOAT_WORD(j,z);
	if (j>0x43000000)				/* if z > 128 */
		return sn*huge*huge;			/* overflow */
	else if (j==0x43000000) {			/* if z == 128 */
		if(p_l+ovt>z-p_h) return sn*huge*huge;	/* overflow */
	}
	else if ((j&0x7fffffff)>0x43160000)		/* z <= -150 */
		return sn*tiny*tiny;			/* underflow */
	else if (j==0xc3160000){			/* z == -150 */
		if(p_l<=z-p_h) return sn*tiny*tiny;		/* underflow */
	}
	/*
	 * compute 2**(p_h+p_l)
	 */
	i = j&0x7fffffff;
	k = (i>>23)-0x7f;
	n = 0;
	if(i>0x3f000000) {		/* if |z| > 0.5, set n = [z+0.5] */
		n = j+(0x00800000>>(k+1));
		k = ((n&0x7fffffff)>>23)-0x7f;	/* new k for n */
		SET_FLOAT_WORD(t,n&~(0x007fffff>>k));
		n = ((n&0x007fffff)|0x00800000)>>(23-k);
		if(j<0) n = -n;
		p_h -= t;
	}
	t = p_l+p_h;
	GET_FLOAT_WORD(is,t);
	SET_FLOAT_WORD(t,is&0xffff8000);
	u = t*lg2_h;
	v = (p_l-(t-p_h))*lg2+t*lg2_l;
	z = u+v;
	w = v-(z-u);
	t = z*z;
	t1 = z - t*(P1+t*(P2+t*(P3+t*(P4+t*P5))));
	r = (z*t1)*__inv(t1-two)-(w+z*w);
	z = one-(r-z);
	GET_FLOAT_WORD(j,z);
	j += (n<<23);
	if((j>>23)<=0) z = __scalbnf(z,n);	/* subnormal output */
	else SET_FLOAT_WORD(z,j);
	return sn*z;
}

/* render a frame of xsz/ysz dimensions into the provided framebuffer */
void render_scanline(int xsz, int ysz, int sl, uint32_t *fb, int samples) {
	int i, s;
	float rcp_samples = __inv((float)samples);

	for(i=0; i<xsz; i++) {
		float r, g, b;
		r = g = b = 0.0f;
			
		for(s=0; s<samples; s++) {
			struct vec3 col = trace(get_primary_ray(i, sl, s), 0);
			r += col.x;
			g += col.y;
			b += col.z;
		}

		r = r * rcp_samples;
		g = g * rcp_samples;
		b = b * rcp_samples;
			
		fb[i] = ((uint32_t)(MIN(r, 1.0f) * 255.0f) & 0xff) << RSHIFT |
					((uint32_t)(MIN(g, 1.0f) * 255.0f) & 0xff) << GSHIFT |
					((uint32_t)(MIN(b, 1.0f) * 255.0f) & 0xff) << BSHIFT | 0xff000000;
	}
}

/* trace a ray throught the scene recursively (the recursion happens through
 * shade() to calculate reflection rays if necessary).
 */
static struct vec3 trace(struct ray ray, int depth) {
	struct vec3 col;
	struct spoint sp, nearest_sp;
	struct sphere *nearest_obj = 0;
	struct sphere *iter = obj_list;

	/* if we reached the recursion limit, bail out */
	if(depth >= MAX_RAY_DEPTH) {
		col.x = col.y = col.z = 0.0f;
		return col;
	}
	
	/* find the nearest intersection ... */
	int j;
	for (j=0; j<nobj; j++) {
		if(ray_sphere(iter, ray, &sp)) {
			if(!nearest_obj || sp.dist < nearest_sp.dist) {
				nearest_obj = iter;
				nearest_sp = sp;
			}
		}
		iter++;
	}

	/* and perform shading calculations as needed by calling shade() */
	if(nearest_obj) {
		col = shade(nearest_obj, &nearest_sp, depth);
	} else {
		col.x = col.y = col.z = 0.0f;
	}

	return col;
}

/* Calculates direct illumination with the phong reflectance model.
 * Also handles reflections by calling trace again, if necessary.
 */
static struct vec3 shade(struct sphere *obj, struct spoint *sp, int depth) {
	int i, j;
	struct vec3 col = {0.0f, 0.0f, 0.0f};

	/* for all lights ... */
	for(i=0; i<lnum; i++) {
		float ispec, idiff;
		struct vec3 ldir;
		struct ray shadow_ray;
		struct sphere *iter = obj_list;
		int in_shadow = 0;

		ldir.x = lights[i].x - sp->pos.x;
		ldir.y = lights[i].y - sp->pos.y;
		ldir.z = lights[i].z - sp->pos.z;

		shadow_ray.orig = sp->pos;
		shadow_ray.dir = ldir;

		/* shoot shadow rays to determine if we have a line of sight with the light */
		for (j=0; j<nobj; j++) {
			if(ray_sphere(iter, shadow_ray, 0)) {
				in_shadow = 1;
				break;
			}
			iter++;
		}

		/* and if we're not in shadow, calculate direct illumination with the phong model. */
		if(!in_shadow) {
			NORMALIZE(ldir);

			idiff = MAX(DOT(sp->normal, ldir), 0.0f);
			ispec = obj->mat.spow > 0.0f ? __powf(MAX(DOT(sp->vref, ldir), 0.0f), obj->mat.spow) : 0.0f;

			col.x += idiff * obj->mat.col.x + ispec;
			col.y += idiff * obj->mat.col.y + ispec;
			col.z += idiff * obj->mat.col.z + ispec;
		}
	}

	/* Also, if the object is reflective, spawn a reflection ray, and call trace()
	 * to calculate the light arriving from the mirror direction.
	 */
	if(obj->mat.refl > 0.0f) {
		struct ray ray;
		struct vec3 rcol;

		ray.orig = sp->pos;
		ray.dir = sp->vref;
		ray.dir.x *= RAY_MAG;
		ray.dir.y *= RAY_MAG;
		ray.dir.z *= RAY_MAG;

		rcol = trace(ray, depth + 1);
		col.x += rcol.x * obj->mat.refl;
		col.y += rcol.y * obj->mat.refl;
		col.z += rcol.z * obj->mat.refl;
	}

	return col;
}

/* calculate reflection vector */
static struct vec3 reflect(struct vec3 v, struct vec3 n) {
	struct vec3 res;
	float dot = v.x * n.x + v.y * n.y + v.z * n.z;
	res.x = -(2.0f * dot * n.x - v.x);
	res.y = -(2.0f * dot * n.y - v.y);
	res.z = -(2.0f * dot * n.z - v.z);
	return res;
}

static struct vec3 cross_product(struct vec3 v1, struct vec3 v2) {
	struct vec3 res;
	res.x = v1.y * v2.z - v1.z * v2.y;
	res.y = v1.z * v2.x - v1.x * v2.z;
	res.z = v1.x * v2.y - v1.y * v2.x;
	return res;
}

/* determine the primary ray corresponding to the specified pixel (x, y) */
static struct ray get_primary_ray(int x, int y, int sample) {
	struct ray ray;
	float m[3][3];
	struct vec3 i, j = {0.0f, 1.0f, 0.0f}, k, dir, orig, foo;

	k.x = cam.targ.x - cam.pos.x;
	k.y = cam.targ.y - cam.pos.y;
	k.z = cam.targ.z - cam.pos.z;
	NORMALIZE(k);

	i = cross_product(j, k);
	j = cross_product(k, i);
	m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
	m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
	m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;
	
	ray.orig.x = ray.orig.y = ray.orig.z = 0.0f;
	ray.dir = get_sample_pos(x, y, sample);
	ray.dir.z = INV_HALF_FOV;
	ray.dir.x *= RAY_MAG;
	ray.dir.y *= RAY_MAG;
	ray.dir.z *= RAY_MAG;
	
	dir.x = ray.dir.x + ray.orig.x;
	dir.y = ray.dir.y + ray.orig.y;
	dir.z = ray.dir.z + ray.orig.z;
	foo.x = dir.x * m[0][0] + dir.y * m[0][1] + dir.z * m[0][2];
	foo.y = dir.x * m[1][0] + dir.y * m[1][1] + dir.z * m[1][2];
	foo.z = dir.x * m[2][0] + dir.y * m[2][1] + dir.z * m[2][2];

	orig.x = ray.orig.x * m[0][0] + ray.orig.y * m[0][1] + ray.orig.z * m[0][2] + cam.pos.x;
	orig.y = ray.orig.x * m[1][0] + ray.orig.y * m[1][1] + ray.orig.z * m[1][2] + cam.pos.y;
	orig.z = ray.orig.x * m[2][0] + ray.orig.y * m[2][1] + ray.orig.z * m[2][2] + cam.pos.z;

	ray.orig = orig;
	ray.dir.x = foo.x + orig.x;
	ray.dir.y = foo.y + orig.y;
	ray.dir.z = foo.z + orig.z;
	
	return ray;
}


static struct vec3 get_sample_pos(int x, int y, int sample) {
	struct vec3 pt;
	static float sf = 0.0f;

	if(sf == 0.0f) {
		sf = 1.5f *__inv((float)xres);
	}

	pt.x = ((float)x * __inv((float)xres)) - 0.5f;
	pt.y = -(((float)y * __inv((float)yres)) - 0.65f) * __inv(aspect);

	if(sample) {
		struct vec3 jt = jitter(x, y, sample);
		pt.x += jt.x * sf;
		pt.y += jt.y * sf * __inv(aspect);
	}
	return pt;
}

/* jitter function taken from Graphics Gems I. */
static struct vec3 jitter(int x, int y, int s) {
	struct vec3 pt;
	pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
	pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
	return pt;
}

/* Calculate ray-sphere intersection, and return {1, 0} to signify hit or no hit.
 * Also the surface point parameters like position, normal, etc are returned through
 * the sp pointer if it is not NULL.
 */
static int ray_sphere(const struct sphere * sph, struct ray ray, struct spoint *sp) {
	float a, b, c, d, sqrt_d, t1, t2;
	
	a = SQ(ray.dir.x) + SQ(ray.dir.y) + SQ(ray.dir.z);
	b = 2.0f * (ray.dir.x * (ray.orig.x - sph->pos.x) +
				ray.dir.y * (ray.orig.y - sph->pos.y) +
				ray.dir.z * (ray.orig.z - sph->pos.z));
	c = SQ(sph->pos.x) + SQ(sph->pos.y) + SQ(sph->pos.z) +
				SQ(ray.orig.x) + SQ(ray.orig.y) + SQ(ray.orig.z) +
				2.0f * (-sph->pos.x * ray.orig.x - sph->pos.y * ray.orig.y - sph->pos.z * ray.orig.z) - SQ(sph->rad);
	
	if((d = SQ(b) - 4.0f * a * c) < 0.0f) return 0;

	sqrt_d = _sqrtf(d);
	float inv2a = 0.5f * __inv(a);
	t1 = (-b + sqrt_d) * inv2a;
	t2 = (-b - sqrt_d) * inv2a;

	if((t1 < ERR_MARGIN && t2 < ERR_MARGIN) || (t1 > 1.0f && t2 > 1.0f)) return 0;

	if(sp) {
		if(t1 < ERR_MARGIN) t1 = t2;
		if(t2 < ERR_MARGIN) t2 = t1;
		sp->dist = t1 < t2 ? t1 : t2;
		
		sp->pos.x = ray.orig.x + ray.dir.x * sp->dist;
		sp->pos.y = ray.orig.y + ray.dir.y * sp->dist;
		sp->pos.z = ray.orig.z + ray.dir.z * sp->dist;
		float invrad = __inv(sph->rad);
		sp->normal.x = (sp->pos.x - sph->pos.x) * invrad;
		sp->normal.y = (sp->pos.y - sph->pos.y) * invrad;
		sp->normal.z = (sp->pos.z - sph->pos.z) * invrad;

		sp->vref = reflect(ray.dir, sp->normal);
		NORMALIZE(sp->vref);
	}
	return 1;
}

static int atomic_fetch_and_increment (unsigned int* global_mutex, int* global_counter)
{
	coprthr_mutex_lock(global_mutex);
	int r = *global_counter;
	int incr = r + 1;
	*global_counter = incr;
	coprthr_mutex_unlock(global_mutex);
	return r;
}

void __entry thread_func(void *tdata)
{
	struct thread_data *td = (struct thread_data*)tdata;

	xres = td->xres;
	yres = td->yres;
	aspect = td->aspect;
	lnum = td->lnum;
	nobj = td->nobj;
	rays_per_pixel = td->rays_per_pixel;

	unsigned int sz = nobj*sizeof(struct sphere) + lnum*sizeof(struct vec3) + sizeof(cam) + NRAN*sizeof(struct vec2) + NRAN*sizeof(int);
	sz = (sz + 7) & 0xfff8; // rounding up to 8 bytes

	void* shared_ptr = coprthr_tls_sbrk(sz);
	int i;
	for (i = 0; i<(sz>>3); i++) {
		((long long*)shared_ptr)[i] = ((long long*)td->shared_ptr)[i];
	}

	obj_list = (struct sphere*)shared_ptr;
	ptrdiff_t offset = nobj*sizeof(struct sphere);
	lights = (struct vec3*)(shared_ptr + offset);
	offset += lnum*sizeof(struct vec3);
	cam = *(struct camera*)(shared_ptr + offset);
	offset += sizeof(struct camera);
	urand = (struct vec2*)(shared_ptr + offset);
	offset += NRAN*sizeof(struct vec2);
	irand = (int*)(shared_ptr + offset);
	offset += NRAN*sizeof(int);

	uint32_t* pixels_ptr = (uint32_t*)(td->shared_ptr + offset);

	unsigned int coreid0 = __corenum_to_ecoreid(0);
	static unsigned int mutex = 0;
	unsigned int* pmutex = (int*)__ega2(coreid0, (void*)&mutex);
	static volatile int counter = 0;
	int* pcounter = (int*)__ega2(coreid0, (void*)&counter);
	int sl = atomic_fetch_and_increment(pmutex, pcounter);

	while (sl < yres) {
		render_scanline(xres, yres, sl, pixels_ptr + sl*xres, rays_per_pixel);
		sl = atomic_fetch_and_increment(pmutex, pcounter);
	}
}
