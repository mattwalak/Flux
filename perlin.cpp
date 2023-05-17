
#include <iostream>
#include "SETTINGS.h"
#include "Util.h"
#include "PerlinNoise.h"
#include "Particle.h"
#include "Color.h"
#include <vector>

#define WIDTH 1920
#define HEIGHT 1080
#define FPS 30.0f
#define START_FRAME 0
#define END_FRAME 30*100
#define NUM_COLORS 6


#define VELOCITY 100.0f
#define TAIL_LENGTH 3.0f
#define STROKE 100.0f

// Animate options
float scale = 5.0f;
float stroke = 100.0f;
float evolution_velocity = 0.25f;
int numParticles = 1000;

// Non-animate options
float step_dt = 2.0f;
int numColors = 6;
std::vector<VEC3> colors;

void addRandomParticle(Tracer * t, float real_time){
	int color_i = (int)rand()%colors.size();
	VEC3 color = colors[color_i];
	VEC2 pos = VEC2(rand()%WIDTH, rand()%HEIGHT);
	t->addParticle(pos, color, real_time, VELOCITY, TAIL_LENGTH, STROKE);
}


int main(){
	// Prepare buffers and solve colors
	float * pixels = new float[WIDTH*HEIGHT*3]();
	solveColors(colors, NUM_COLORS);
	float real_time;

	// Prepare simulation
	Tracer * t = new Tracer(WIDTH, HEIGHT, scale, evolution_velocity);

	for(int frame_num = START_FRAME; frame_num < END_FRAME; frame_num++){
		std::cout << "rendering frame #" << frame_num << std::endl;
		real_time = frame_num*(1/FPS);

		// ADD PARTICLES AND UPDATE ANIMATON DATA HERE vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		for(int i = 0; i < 10; i++)
			addRandomParticle(t, real_time);
		
		// ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

		t->render(real_time, pixels);
		char filename[100];
		sprintf(filename, "oven/frame_%04i.ppm", frame_num);
		writePPM(filename, WIDTH, HEIGHT, pixels);
	}
}
