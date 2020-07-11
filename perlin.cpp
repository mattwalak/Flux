
#include <iostream>
#include "SETTINGS.h"
#include "Util.h"
#include "PerlinNoise.h"
#include "Particle.h"
#include "hsvtorgb.h"
#include <vector>

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

void solveColors(std::vector<VEC3>& colors, int numColors){
	float component_sum = 0.0f;
	for(int i = 0; i < numColors; i++){
		float thisColor[3];
		float h = 360.0f*(float)i/numColors;
		HSVtoRGB(h, 1.0f, 1.0f, thisColor);
		component_sum += thisColor[0];
		colors.push_back(VEC3(thisColor[0], thisColor[1], thisColor[2]));
	}
	for(int i = 0; i < numColors*3; i++){
		colors[i] /= component_sum;
	}
}

int main(){
	std::vector<VEC3> colors;
	solveColors(colors, 6);
	Tracer t(width, height, 10000, scale, colors);
	for(int i = 0; i < 1000; i++){
		if(i%100 == 0){
			int count = t.particleCount();
			std::cout << i << ") num particles = " << count << "\n";
			if(count == 0)
				break;
		}
		t.updateVelocity();
		t.step();
	}
	float * pixels = new float[width*height*3]();
	t.toImage(pixels);
	writePPM("line.ppm", width, height, pixels);

	/*
    pNoise = new PerlinNoise(2048);
    float * pixels = new float[width*height*3];
    genTexture(pixels);
    writePPM("noise.ppm", width, height, pixels);*/
}
