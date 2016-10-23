#ifndef _common_h_
#define _common_h_

struct vec2 {
	float x, y;
} __attribute__((packed));

struct vec3 {
	float x, y, z;
} __attribute__((packed));

struct ray {
	struct vec3 orig, dir;
} __attribute__((packed));

struct material {
	struct vec3 col;	/* color */
	float spow;		/* specular power */
	float refl;		/* reflection intensity */
} __attribute__((packed));

struct sphere {
	struct vec3 pos;
	float rad;
	struct material mat;
} __attribute__((packed));

struct spoint {
	struct vec3 pos, normal, vref;	/* position, normal and view reflection */
	float dist;		/* parametric distance of intersection along the ray */
} __attribute__((packed));

struct camera {
	struct vec3 pos, targ;
	float fov;
} __attribute__((packed));

struct thread_data {
	int xres;
	int yres;
	float aspect;
	int nobj;
	int lnum;
	int rays_per_pixel;
	void* shared_ptr;
} __attribute__((packed));

/* global state */
int xres = 800;
int yres = 600;
float aspect;
struct sphere *obj_list;
int nobj;
struct vec3 *lights;
int lnum;
int rays_per_pixel;
struct camera cam;

#define NRAN 64
#define MASK (NRAN - 1)
struct vec2* urand;
int* irand;

/* bit-shift ammount for packing each color into a 32bit uint */
#define RSHIFT 16
#define GSHIFT 8
#define BSHIFT 0

#endif // _common_h_
