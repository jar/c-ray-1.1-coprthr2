/* c-ray-mt - a simple multithreaded raytracing filter.
 * Copyright (C) 2006 John Tsiombikas <nuclear@siggraph.org>
 *
 * You are free to use, modify and redistribute this program under the
 * terms of the GNU General Public License v2 or (at your option) later.
 * see "http://www.gnu.org/licenses/gpl.txt" for details.
 * ---------------------------------------------------------------------
 * Usage:
 *   compile:  just type make
 *              (add any arch-specific optimizations for your compiler in CFLAGS first)
 *       run:  cat scene | ./c-ray-mt [-t num-threads] >foo.ppm
 *              (on broken systems such as windows try: c-ray-mt -i scene -o foo.ppm)
 *     enjoy:  display foo.ppm
 *              (with imagemagick, or use your favorite image viewer)
 * ---------------------------------------------------------------------
 * Scene file format:
 *   # sphere (many)
 *   s  x y z  rad   r g b   shininess   reflectivity
 *   # light (many)
 *   l  x y z
 *   # camera (one)
 *   c  x y z  fov   tx ty tz
 * ---------------------------------------------------------------------
 */

/*
 * Modified to support the Adapteva Epiphany architecure using the COPRTHR-2
 * interface.
 * Copyright (C) 2016 James A. Ross (james.a.ross@gmail.com)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#ifdef USESDL
#include <pthread.h>
#include <SDL.h>
#endif
#include "coprthr.h"
#include "coprthr_cc.h"
#include "coprthr_thread.h"
#include "common.h"

#define XSTR(s) STR(s)
#define STR(s) #s
#define DIE(...) do{fprintf(stderr,__VA_ARGS__); exit(EXIT_FAILURE);}while(0)

#define VER_MAJOR	1
#define VER_MINOR	1
#define VER_STR "c-ray-mt v"XSTR(VER_MAJOR)"."XSTR(VER_MINOR)" ... for Epiphany!\n"

void load_scene(FILE *fp);
unsigned long get_msec(void);

const char *usage = {
	"Usage: c-ray-mt [options]\n"
	"  Reads a scene file from stdin, writes the image to stdout, and stats to stderr.\n\n"
	"Options:\n"
	"  -t <num>   how many threads to use (default: 1)\n"
	"  -s WxH     where W is the width and H the height of the image\n"
	"  -r <rays>  shoot <rays> rays per pixel (antialiasing)\n"
	"  -i <file>  read from <file> instead of stdin\n"
	"  -o <file>  write to <file> instead of stdout\n"
	"  -h         this help screen\n\n"
};

#ifdef USESDL
void *SDL_thread(void *tdata)
{
	struct thread_data *td = (struct thread_data*)tdata;
	int xres = td->xres;
	int yres = td->yres;
	int lnum = td->lnum;
	void* shared_ptr = (void*)td->shared_ptr;
	ptrdiff_t offset = nobj*sizeof(struct sphere) + lnum*sizeof(struct vec3) + sizeof(cam) + NRAN*sizeof(struct vec2) + NRAN*sizeof(int);

	uint32_t* pixels = (uint32_t*)(shared_ptr+offset);

	SDL_Window* window;
	SDL_Renderer* renderer;
	SDL_Surface* surface;
	SDL_Texture* texture;
	SDL_Event event;
	SDL_bool done = SDL_FALSE;

	SDL_LogSetPriority(SDL_LOG_CATEGORY_APPLICATION, SDL_LOG_PRIORITY_INFO);

	if (SDL_Init(SDL_INIT_VIDEO) < 0) DIE("Couldn't initialize SDL: %s\n", SDL_GetError());

	window = SDL_CreateWindow(VER_STR, SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, xres, yres, SDL_WINDOW_RESIZABLE);
	if (!window) DIE("Couldn't create window: %s\n", SDL_GetError());

	renderer = SDL_CreateRenderer(window, -1, 0);
	if (!renderer) DIE("Couldn't create renderer: %s\n", SDL_GetError());

	unsigned int rmask = 0xff << RSHIFT;
	unsigned int gmask = 0xff << GSHIFT;
	unsigned int bmask = 0xff << BSHIFT;
	unsigned int amask = 0xff << 24;
	int pitch = xres*sizeof(unsigned int);

	while (!done) {
	surface = SDL_CreateRGBSurfaceFrom(pixels, xres, yres, 32, pitch, rmask, gmask, bmask, amask);
	if (!surface) fprintf(stderr, "Couldn't create surface: %s\n", SDL_GetError());

	texture = SDL_CreateTextureFromSurface(renderer, surface);
	if (!texture) DIE("Couldn't create texture: %s\n", SDL_GetError());

		while (SDL_PollEvent(&event)) {
			switch (event.type) {
				case SDL_KEYDOWN:
					if (event.key.keysym.sym == SDLK_ESCAPE)
						done = SDL_TRUE;
					break;
				case SDL_QUIT:
					done = SDL_TRUE;
					break;
			}
		}

		SDL_RenderClear(renderer);
		SDL_RenderCopy(renderer, texture, NULL, NULL);
		SDL_RenderPresent(renderer);
		SDL_Delay(200);

	SDL_DestroyTexture(texture);
	SDL_FreeSurface(surface);
	}
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	pthread_exit(NULL);
}
#endif // USESDL

int main(int argc, char **argv)
{
	int i;
	unsigned long rend_time, start_time;
	uint32_t *pixels;
	float sl;
	int rays_per_pixel = 1;
	int thread_num = 1;
	FILE *infile = stdin, *outfile = stdout;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-' && argv[i][2] == 0) {
			char *sep;
			switch(argv[i][1]) {
			case 't':
				if(!isdigit(argv[++i][0]))
					DIE("-t mus be followed by the number of worker threads to spawn\n");
				thread_num = atoi(argv[i]);
				if(!thread_num)
					DIE("invalid number of threads specified: %d\n", thread_num);
				break;
					
			case 's':
				if(!isdigit(argv[++i][0]) || !(sep = strchr(argv[i], 'x')) || !isdigit(*(sep + 1)))
					DIE("-s must be followed by something like \"640x480\"\n");
				xres = atoi(argv[i]);
				yres = atoi(sep + 1);
				aspect = (float)xres / (float)yres;
				break;

			case 'i':
				if(!(infile = fopen(argv[++i], "rb")))
					DIE("failed to open input file %s: %s\n", argv[i], strerror(errno));
				break;

			case 'o':
				if(!(outfile = fopen(argv[++i], "wb")))
					DIE("failed to open output file %s: %s\n", argv[i], strerror(errno));
				break;

			case 'r':
				if(!isdigit(argv[++i][0]))
					DIE("-r must be followed by a number (rays per pixel)\n");
				rays_per_pixel = atoi(argv[i]);
				break;

			case 'h':
				fputs(usage, stdout);
				return 0;
				
			default:
				DIE("unrecognized argument: %s\n%s", argv[i], usage);
			}
		}
		else
			DIE("unrecognized argument: %s\n%s", argv[i], usage);
	}

	aspect = (float)xres/(float)yres;

	load_scene(infile);

	/* initialize the random number tables for the jitter */
	urand = (struct vec2*)malloc(NRAN*sizeof(struct vec2));
	irand = (int*)malloc(NRAN*sizeof(int));
	for(i=0; i<NRAN; i++) urand[i].x = (float)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) urand[i].y = (float)rand() / RAND_MAX - 0.5;
	for(i=0; i<NRAN; i++) irand[i] = (int)(NRAN * ((float)rand() / RAND_MAX));

	if(thread_num > yres) {
		fprintf(stderr, "more threads than scanlines specified, reducing number of threads to %d\n", yres);
		thread_num = yres;
	}

	int dd = coprthr_dopen(COPRTHR_DEVICE_E32,COPRTHR_O_THREAD);
	if (dd<0) DIE("device open failed\n");

	coprthr_program_t prg = coprthr_cc_read_bin("device.e32", 0);
	if (!prg) DIE("Unable to read program\n");

	coprthr_sym_t thr = coprthr_getsym(prg,"thread_func");
	if (!thr) DIE("Unable to get thread function symbol\n");

	size_t sz = nobj*sizeof(struct sphere) + lnum*sizeof(struct vec3) + sizeof(cam) + NRAN*sizeof(struct vec2) + NRAN*sizeof(int) + xres*yres*sizeof(*pixels);
	ptrdiff_t offset = 0;
	coprthr_mem_t shared_mem = coprthr_dmalloc(dd, sz, 0);
	void* shared_ptr = coprthr_memptr(shared_mem, 0);
	bzero(shared_ptr, sz);

	memcpy(shared_ptr+offset, (void*)obj_list, nobj*sizeof(struct sphere));
	offset += nobj*sizeof(struct sphere);
	memcpy(shared_ptr+offset, (void*)lights, lnum*sizeof(struct vec3));
	offset += lnum*sizeof(struct vec3);
	memcpy(shared_ptr+offset, (void*)&cam, sizeof(cam));
	offset += sizeof(cam);
	memcpy(shared_ptr+offset, (void*)urand, NRAN*sizeof(struct vec2));
	offset += NRAN*sizeof(struct vec2);
	memcpy(shared_ptr+offset, (void*)irand, NRAN*sizeof(int));
	offset += NRAN*sizeof(int);

	coprthr_td_t td;
	coprthr_attr_t attr;
	coprthr_attr_init(&attr);
	coprthr_attr_setdetachstate(&attr,COPRTHR_CREATE_JOINABLE);
	coprthr_attr_setdevice(&attr,dd);
	void* status;
	
	fprintf(stderr, VER_STR);
	
	struct thread_data args = {
		.xres = xres,
		.yres = yres,
		.aspect = aspect,
		.nobj = nobj,
		.lnum = lnum,
		.rays_per_pixel = rays_per_pixel,
		.shared_ptr = shared_ptr
	};

	fprintf(stderr, "threads: %d, args: xres = %d, yres = %d, aspect = %f, nobj = %d, lnum = %d, rays_per_pixel = %d, shared_ptr (%d bytes) = %p\n", thread_num, args.xres, args.yres, args.aspect, args.nobj, args.lnum, args.rays_per_pixel, sz, args.shared_ptr);

	coprthr_mem_t argmem = coprthr_dmalloc(dd,sizeof(args),0);
	coprthr_dwrite(dd,argmem,0,&args,sizeof(args),COPRTHR_E_NOW);

#ifdef USESDL
	pthread_t sdl_thread;
	int rc = pthread_create(&sdl_thread, NULL, SDL_thread, (void *)&args);
	if (rc) DIE("ERROR; return code from pthread_create() is %d\n", rc);
#endif

	start_time = get_msec();

	coprthr_ncreate(thread_num, &td, &attr, thr, (void*)&argmem);
	coprthr_attr_destroy(&attr);
	coprthr_join(td,&status);

	rend_time = get_msec() - start_time;
	
	pixels = (uint32_t*)(shared_ptr+offset);

	/* output statistics to stderr */
	fprintf(stderr, "Rendering took: %lu seconds (%lu milliseconds)\n", rend_time / 1000, rend_time);

	/* output the image */
	fprintf(outfile, "P6\n%d %d\n255\n", xres, yres);
	for(i=0; i<xres * yres; i++) {
		fputc((pixels[i] >> RSHIFT) & 0xff, outfile);
		fputc((pixels[i] >> GSHIFT) & 0xff, outfile);
		fputc((pixels[i] >> BSHIFT) & 0xff, outfile);
	}
	fflush(outfile);

#ifdef USESDL
	rc = pthread_join(sdl_thread, NULL);
#endif

	coprthr_dfree(dd,argmem);
	coprthr_dfree(dd,shared_mem);

	if(infile != stdin) fclose(infile);
	if(outfile != stdout) fclose(outfile);
	return 0;
}

/* Load the scene from an extremely simple scene description file */
#define DELIM	" \t\n"
void load_scene(FILE *fp) {
	char line[256], *ptr, type;

	while((ptr = fgets(line, 256, fp))) {
		int i;
		struct vec3 pos, col;
		float rad, spow, refl;
		
		while(*ptr == ' ' || *ptr == '\t') ptr++;
		if(*ptr == '#' || *ptr == '\n') continue;

		if(!(ptr = strtok(line, DELIM))) continue;
		type = *ptr;
		
		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((float*)&pos.x + i) = atof(ptr);
		}

		if(type == 'l') {
			lnum++;
			lights = (struct vec3*)realloc(lights, lnum*sizeof(*lights));
			lights[lnum-1] = pos;
			continue;
		}

		if(!(ptr = strtok(0, DELIM))) continue;
		rad = atof(ptr);

		for(i=0; i<3; i++) {
			if(!(ptr = strtok(0, DELIM))) break;
			*((float*)&col.x + i) = atof(ptr);
		}

		if(type == 'c') {
			cam.pos = pos;
			cam.targ = col;
			cam.fov = rad;
			continue;
		}

		if(!(ptr = strtok(0, DELIM))) continue;
		spow = atof(ptr);

		if(!(ptr = strtok(0, DELIM))) continue;
		refl = atof(ptr);

		if(type == 's') {
			nobj++;
			obj_list = (struct sphere*)realloc(obj_list, nobj*sizeof(*obj_list));
			struct sphere* sph = &obj_list[nobj - 1];
			sph->pos = pos;
			sph->rad = rad;
			sph->mat.col = col;
			sph->mat.spow = spow;
			sph->mat.refl = refl;
		} else {
			fprintf(stderr, "unknown type: %c\n", type);
		}
	}
}

unsigned long get_msec(void) {
	static struct timeval timeval, first_timeval;
	
	gettimeofday(&timeval, 0);
	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}
