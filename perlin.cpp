
#include <iostream>
#include "SETTINGS.h"
#include "Util.h"
#include "PerlinNoise.h"
#include "Particle.h"
#include "Color.h"
#include <vector>

#define WIDTH 640
#define HEIGHT 360
#define FPS 30.0f
#define START_FRAME 0
#define END_FRAME 300
#define NUM_COLORS 6

// Animate options
float scale = 2.0f;
int stroke = 2;
float evolution = 0.0f;
int numParticles = 1000;

// Non-animate options
float step_dt = 2.0f;
int numColors = 6;


int main(){
	// Prepare buffers and solve colors
	float * pixels = new float[WIDTH*HEIGHT*3]();
	float real_time;
	std::vector<VEC3> colors;
	solveColors(colors, NUM_COLORS);

	// Prepare simulation
	Tracer * t = new Tracer(WIDTH, HEIGHT, scale, colors);

	for(int frame_num = START_FRAME; frame_num < END_FRAME; frame_num++){
		std::cout << "rendering frame #" << frame_num << std::endl;
		real_time = frame_num*(1/FPS);
		
		// --------------------------------------------------------------

		t->setRealTime(real_time);

		// Add particles (If applicable) and update anim information
		if(frame_num == 2){
			t->addParticle(real_time, 100.0f, 1.0f);
		}
		
		


		t->resetBuffers();
		t->simulateFrame(0.0f);
		t->toImage(pixels);

		// --------------------------------------------------------------

		char filename[100];
		sprintf(filename, "oven/frame_%04i.ppm", frame_num);
		writePPM(filename, WIDTH, HEIGHT, pixels);
	}
}
