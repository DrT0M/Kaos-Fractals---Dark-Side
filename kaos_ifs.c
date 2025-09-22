//	Kaos IFS (Iterated Function Systems) fractals

//Compile:
//	cc kaos_ifs.c -o ~/kaos_ifs -lX11 -lm

//Run with the default parameters [9000 1000 0.6]:
//	~/kaos_ifs

//Run with specific parameters:
//	~/kaos_ifs 20000 700 0.55

//  tab stop @ 8
/*****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
//	#include "vroot.h"	//XScreenSaver

/*****************************************************************************/
	Display*display;
	Window  window;
	Visual *visual;
	GC 	gc;
	XImage *ximage;
	unsigned long*_rImage;
	unsigned long**rImage;
	unsigned long*_gImage;
	unsigned long**gImage;
	unsigned long*_bImage;
	unsigned long**bImage;
	double*A11;
	double*A12;
	double*A21;
	double*A22;
	double*B1;
	double*B2;
	double*RG;
	double*GB;
	double x, xmin, xmax, xmid;	int xs;
	double y, ymin, ymax, ymid;	int ys; 
	double maxmid, xOld, tscale, r;
	double seed, oldseed;
	int Nn, Nmax, opsmax, jop, jopmax;
	int i, j, jOld1, jOld2;
	int glowing, FadeOut, EndProgram, MaxScreens, FullScreen, Clip;
	unsigned long MyR, MyG, MyB;
	unsigned long MyRmin, MyGmin, MyBmin;
	unsigned long MyRmax, MyGmax, MyBmax;
	int argv1, argv2; float argv3;

#define	TRUE	(1==1)
#define	FALSE	(1==0)

/*****************************************************************************/
//	splitmix64.c
/* This is a fixed-increment version of Java 8's SplittableRandom generator
   See http://dx.doi.org/10.1145/2714064.2660195 and 
   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

   It is a very fast generator passing BigCrush, and it can be useful if
   for some reason you absolutely want 64 bits of state. */
/*
	Initialization
	It is the recommendation of the authors of the xoshiro paper 
	to initialize the state of the generators using a generator 
	which is radically different from the initialized generators, 
	as well as one which will never give the "all-zero" state; 
	for shift-register generators, this state is impossible to escape from.[14][15] 
	The authors specifically recommend using the SplitMix64 generator, from a 64-bit seed, as follows:

#include <stdint.h>
 */
struct splitmix64_state 
{
	uint64_t s;
};

uint64_t splitmix64(struct splitmix64_state *state) 
{
	uint64_t result = state->s;

	state->s = result + 0x9E3779B97f4A7C15;
	result = (result ^ (result >> 30)) * 0xBF58476D1CE4E5B9;
	result = (result ^ (result >> 27)) * 0x94D049BB133111EB;
	return result ^ (result >> 31);
}

/*****************************************************************************/
/*
//	64bits_xoshiro256plus.c
//	This is xoshiro256+ 1.0, our best and fastest generator for floating-point

	xoshiro256+ is approximately 15% faster than xoshiro256**, 
	but the lowest three bits have low linear complexity; 
	therefore, it should be used only for floating point results by extracting the upper 53 bits.

	The state must be seeded so that it is not everywhere zero. If you have
	a 64-bit seed, we suggest to seed a splitmix64 generator and use its
	output to fill s. 

#include <stdint.h>
 */
static inline uint64_t rotl(const uint64_t x, int k) 
{
	return (x << k) | (x >> (64 - k));
}

static uint64_t state_xoshiro256p[4];

uint64_t next_xoshiro256p(void) 
{
	uint64_t*s = state_xoshiro256p;

//	const uint64_t result = rotl(s[0] + s[3], 23) + s[0];	//xoshiro256plusplus
	const uint64_t result = s[0] + s[3];	//xoshiro256plus

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

void init_xoshiro256p(uint64_t seed) 
{
	state_xoshiro256p[0] = 0;
	state_xoshiro256p[1] = 0;
	state_xoshiro256p[2] = 0;
	state_xoshiro256p[3] = 0;

	struct splitmix64_state smstate = {seed};

	uint64_t tmp = splitmix64(&smstate);
	state_xoshiro256p[0] = (uint32_t)tmp;
	state_xoshiro256p[1] = (uint32_t)(tmp >> 32);

	tmp = splitmix64(&smstate);
	state_xoshiro256p[2] = (uint32_t)tmp;
	state_xoshiro256p[3] = (uint32_t)(tmp >> 32);
}

/*****************************************************************************/
//	Kaos IFS (Iterated Function Systems) fractals is licensed under
//	MIT License
//
//	Copyright (c) 2025 DrT0M
//
//	Permission is hereby granted, free of charge, to any person obtaining a copy
//	of this software and associated documentation files (the "Software"), to deal
//	in the Software without restriction, including without limitation the rights
//	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//	copies of the Software, and to permit persons to whom the Software is
//	furnished to do so, subject to the following conditions:
//	
//	The above copyright notice and this permission notice shall be included in all
//	copies or substantial portions of the Software.
//
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//	SOFTWARE.
//

/*****************************************************************************/
double uniform(double lowerbound, double upperbound)
{
	return lowerbound + (upperbound-lowerbound)*(((next_xoshiro256p(),next_xoshiro256p()) >> 11) * 0x1.0p-53);
}

/*****************************************************************************/
typedef union
{
	uint64_t i;
	double d;
}
U;

void Randomize( int now, U seed )
{
	if	( ! now )
	{
		init_xoshiro256p( seed.i );
	}
	else
	{
		struct timeval tv;
		(void) gettimeofday(&tv, NULL);
		init_xoshiro256p( tv.tv_sec*1000000 + tv.tv_usec );
	}
}

/*****************************************************************************/
static double lastRnd = 0.;

double Rnd( int action, U seed )
{
//	action < 0	The same number every time, using seed.
//	action > 0	The next number in the pseudo-random sequence.
//	action == 0	The most recently generated number.

	if	( action < 0 )	
		init_xoshiro256p( seed.i );

	if	( action != 0 )	lastRnd = uniform( 0.0, 1.0 ); 

	return	lastRnd;
}

/*****************************************************************************/
void PSet(Display*display, Window window, GC gc, int h, int v, unsigned long rgb)
{
	XSetForeground(	display		,gc ,rgb	);
	XDrawPoint(	display	,window	,gc	,h ,v	);
}

/*****************************************************************************/
unsigned long setRGB(double R, double G, double B)
{
	return	(((int)(R*255)%256)<<16)+
		(((int)(G*255)%256)<<8 )+
		(((int)(B*255)%256)    );
}

/*****************************************************************************/
unsigned long RGB(unsigned long R, unsigned long G, unsigned long B)
{
	return	(((int)(R)%256)<<16)+
		(((int)(G)%256)<<8 )+
		(((int)(B)%256)    );
}

/*****************************************************************************/
void Initialize( double seedling )
{
	for( i = 0 ; i < 1000 ; i++ )	(void)	Rnd(1,(U)0.);

	r = Rnd(-1, (U)seedling);
	Randomize( FALSE, (U)r );

	Nn = (int)(Rnd(1,(U)0.) * (Nmax - 4)) + 5; 
	if (Nn < 2)
		Nn = 2;

	A11 = realloc((A11),0);		A11 = calloc( Nn, sizeof( double ) );
	A12 = realloc((A12),0);		A12 = calloc( Nn, sizeof( double ) );
	A21 = realloc((A21),0);		A21 = calloc( Nn, sizeof( double ) );
	A22 = realloc((A22),0);		A22 = calloc( Nn, sizeof( double ) );
	B1  = realloc((B1 ),0);		B1  = calloc( Nn, sizeof( double ) );
	B2  = realloc((B2 ),0);		B2  = calloc( Nn, sizeof( double ) );
	RG  = realloc((RG ),0);		RG  = calloc( Nn, sizeof( double ) );
	GB  = realloc((GB ),0);		GB  = calloc( Nn, sizeof( double ) );

	int i;
	for( i = 0 ; i < Nn ; i++ )
	{
		A11[i] = tscale;
		A22[i] = tscale;
		A12[i] = (Rnd(1,(U)0.) * 2 - 1) * tscale;
		A21[i] = (Rnd(1,(U)0.) * 2 - 1) * tscale;
		B1[i] = (Rnd(1,(U)0.) * 2 - 1) * tscale;
		B2[i] = (Rnd(1,(U)0.) * 2 - 1) * tscale;
		RG[i] = Rnd(1,(U)0.);
		r = Rnd(1,(U)0.);
		if (RG[i] < r )
		{
			GB[i] = r;
		}
		else
		{
			GB[i] = RG[i];
			RG[i] = r;
		}
	}
	x = 0;
	y = 0;
	for( i = 0 ; i < 100 ; i++ )
	{
		jOld1 = j;
		j = floor(Rnd(1,(U)0.) * Nn) + 1;
		xOld = x;
		x = A11[j] * xOld + A12[j] * y + B1[j];
		y = A21[j] * xOld + A22[j] * y + B2[j];
	}
	jop = 0;

}

/*****************************************************************************/
#define	fading_TRUE	(1==1)
#define	fading_FALSE	(1==0)

#define Fade()	Glowing(fading_TRUE);
#define Glow()	Glowing(fading_FALSE);

/*****************************************************************************/
void Glowing( int fading )
{
	jOld2 = jOld1;
	jOld1 = j;
	r = Rnd(1,(U)0.);
	j = floor(r * Nn);
	xOld = x;
	x = A11[j] * xOld + A12[j] * y + B1[j];
	y = A21[j] * xOld + A22[j] * y + B2[j];
	xs = floor(xmid + x * maxmid);
	ys = floor(ymid + y * maxmid);

	if( (xs != floor(xmid)) || (ys != floor(ymid)) ) 
	if( (xs >= 0) && (ys >= 0) && (xs < xmax) && (ys < ymax) )
	{
		MyR = rImage[ xs ][ ys ];
		MyG = gImage[ xs ][ ys ];
		MyB = bImage[ xs ][ ys ];

		//FadeOut -1:Reverse, 0:Blackout, 1:Clear Screen

		if	( r < RG[ (floor(r * 3) == 0. ? j : (floor(r * 3) == 1. ? jOld1 : jOld2)) ] )
		{
			if( fading )
			{
				MyR = (MyR <= MyRmin ? 0UL : (MyRmax - (MyRmax - MyR) * 2UL));
				if( FadeOut == 0 )
					MyR = 0UL;
			}
			else
			{
				MyR = (MyR <= MyRmin ? MyRmin : (MyRmax - (MyRmax - MyR) / 2UL));
			}
		}
		else if	( r < GB[ (floor(r * 3) == 0. ? j : (floor(r * 3) == 1. ? jOld1 : jOld2)) ] )
		{
			if( fading )
			{
				MyG = (MyG <= MyGmin ? 0UL : (MyGmax - (MyGmax - MyG) * 2UL));
				if( FadeOut == 0 )
					MyG = 0UL;
			}
			else
			{
				MyG = (MyG <= MyGmin ? MyGmin : (MyGmax - (MyGmax - MyG) / 2UL));
			}
		}
		else
		{
			if( fading )
			{
				MyB = (MyB <= MyBmin ? 0UL : (MyBmax - (MyBmax - MyB) * 2UL));
				if( FadeOut == 0 )
					MyB = 0UL;
			}
			else
			{
				MyB = (MyB <= MyBmin ? MyBmin : (MyBmax - (MyBmax - MyB) / 2UL));
			}
		}

		rImage[ xs ][ ys ] = MyR;
		gImage[ xs ][ ys ] = MyG;
		bImage[ xs ][ ys ] = MyB;

		PSet(display, window, gc, xs, ys, RGB(MyR,MyG,MyB));
	}
}

/*****************************************************************************/
void Timer1_Timer()
{
	if( glowing )
		usleep (argv1);
	else
		usleep (argv2);

	int i;
	for( i = 0 ; i < opsmax ; i++ )
	{
		jop = jop + 1;
		if( glowing )
		{
			Glow();
			if( jop > xmax * ymax * MaxScreens || MyR == MyRmax || MyG == MyGmax || MyB == MyBmax )
			{
				oldseed = seed;
				jopmax = jop; 
				if( FadeOut == 1 ) jopmax = 0;	//FadeOut -1:Reverse, 0:Blackout, 1:Clear Screen
				glowing = FALSE;
				Initialize(seed);
				break;
			}
		}
		else
		{
			Fade();
			if( jop > jopmax )
			{
				glowing = TRUE;
				XClearWindow(display, window);
				oldseed = seed;
				seed = Rnd(1,(U)0.);
				Initialize(seed);
				break;
			}
		}
	}
}

/*****************************************************************************/
int main(int argc, char*argv[])
{
	argv1 = (argc > 1 ? atoi(argv[1]) : 9000);
	argv2 = (argc > 2 ? atoi(argv[2]) : 1000);
	argv3 = (argc > 3 ? atof(argv[3]) : 0.60);

	/*****************************************************************************/
	display = XOpenDisplay(getenv("DISPLAY"));
	window  = DefaultRootWindow(display);
	visual  = DefaultVisual(display, 0);
	if(visual->class != TrueColor)
	{
		fprintf(stderr, "Cannot handle non-true color visual ...\n");
		exit(1);
	}
	XWindowAttributes window_attr;
	XGetWindowAttributes(display, window, &window_attr);
	int width  = window_attr.width ;
	int height = window_attr.height;

	/*****************************************************************************/
#/*[*/if !defined(_VROOT_H_)

	display = XOpenDisplay(NULL);
	window  = XCreateSimpleWindow(display, RootWindow(display, 0), 0, 0, width, height, 1, 0, 0);
	visual  = DefaultVisual(display, 0);
	if(visual->class != TrueColor)
	{
		fprintf(stderr, "Cannot handle non-true color visual ...\n");
		exit(1);
	}

#/*]*/endif

	/*****************************************************************************/
	gc = XCreateGC(display, window, 0, NULL);
	XSetForeground(display, gc, WhitePixelOfScreen(DefaultScreenOfDisplay(display)));

	char*image32 = (char*) calloc(width*height*4, sizeof(char));
	ximage = XCreateImage(display, visual, 24, ZPixmap, 0, image32, width, height, 32, 0);

	/*****************************************************************************/
#/*[*/if !defined(_VROOT_H_)

	for( int k=0 ; k<argc-1 ; k++ )		argv[k][ strlen(argv[k]) ] = ' ';
	XStoreName(display, window, argv[0]);
	XSelectInput(display, window, ButtonPressMask | ExposureMask | StructureNotifyMask);
	XMapWindow(display, window);
	for(;;)
	{
		XEvent ev;
		XNextEvent(display, &ev);
		if(ev.type == MapNotify)	break;
	}
	//Wait for the MapNotify event
	//which means that the window has appeared on the screen.
	//initial black screen
	XPutImage(display, window, DefaultGC(display, 0), ximage, 0, 0, 0, 0, width, height);
	XFlush(display);

#/*]*/endif

	/*****************************************************************************/
/*	unsigned long*	*/	_rImage = (unsigned long*) calloc( width * height, sizeof( unsigned long ) );
/*	unsigned long**	*/	 rImage = (unsigned long**)malloc( width *         sizeof( unsigned long*) );
	for( int k = 0 ; k < width; k++ )	rImage[ k ] = _rImage + k * height;

/*	unsigned long*	*/	_gImage = (unsigned long*) calloc( width * height, sizeof( unsigned long ) );
/*	unsigned long**	*/	 gImage = (unsigned long**)malloc( width *         sizeof( unsigned long*) );
	for( int k = 0 ; k < width; k++ )	gImage[ k ] = _gImage + k * height;

/*	unsigned long*	*/	_bImage = (unsigned long*) calloc( width * height, sizeof( unsigned long ) );
/*	unsigned long**	*/	 bImage = (unsigned long**)malloc( width *         sizeof( unsigned long*) );
	for( int k = 0 ; k < width; k++ )	bImage[ k ] = _bImage + k * height;

	/*****************************************************************************/
	FullScreen = 1;
	EndProgram = FALSE;
	MyRmin=97;
	MyRmax=222;
	MyGmin=96;
	MyGmax=223;
	MyBmin=128;
	MyBmax=254;
	Clip = (0 < 1);
	FadeOut=-1;	//FadeOut -1:Reverse, 0:BlackOut, 1:ClearScreen
	tscale = argv3;
	MaxScreens = 1;
	opsmax = 1000;
	Nmax = 18;

	xmin = 0; xmax = width ; xmid = xmax / 2;
	ymin = 0; ymax = height; ymid = ymax / 2;
	maxmid = (ymid > xmid ? ymid : xmid);

	XClearWindow(display, window);

	Randomize( TRUE, (U)0. );
	glowing = TRUE;
	oldseed = 0;
	seed = Rnd(1,(U)0.);
	Initialize(seed);


	/*****************************************************************************/
	while (1)
	{
		Timer1_Timer();
		XFlush(display);
	}
	XCloseDisplay(display);
}
