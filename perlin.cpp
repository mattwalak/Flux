
#include <iostream>
#include "SETTINGS.h"
#include "Util.h"
#include "PerlinNoise.h"
#include "Particle.h"

int width = 1000;
int height = 1000;
float scale = 2.0f;
PerlinNoise * pNoise;
Particle * particles;

void genTexture(float * pixels){
	for(int x = 0; x < width; x++){
		for(int y = 0; y < height; y++){
			float u = (float)x/width;
			float v = (float)y/height;
			float sample = pNoise->noise(u*scale, v*scale);
			sample += 1.0f;
			sample *= 0.5;
			sample = clamp(sample);

			// std::cout << sample << std::endl;

			pixels[3*(y*width + x) + 0] = sample * 255.0f;
			pixels[3*(y*width + x) + 1] = sample * 255.0f;
			pixels[3*(y*width + x) + 2] = sample * 255.0f;
		}
	}
}

int main(){
	Tracer t(width, height, 1000000, scale);
	for(int i = 0; i < 100000; i++){
		if(i%100 == 0){
			int count = t.particleCount();
			std::cout << i << ") num particles = " << count << "\n";
			if(count == 0)
				break;
		}
		t.updateVelocity();
		t.step();
	}
	float * pixels = new float[width*height*3];
	t.toImage(pixels);

	writePPM("line.ppm", width, height, pixels);

	/*
    pNoise = new PerlinNoise(2048);
    float * pixels = new float[width*height*3];
    genTexture(pixels);
    writePPM("noise.ppm", width, height, pixels);*/
}
