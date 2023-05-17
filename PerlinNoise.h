#ifndef PerlinNoise_h
#define PerlinNoise_h

#include "SETTINGS.h"
#include "Util.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

#define LERP(start, end, t) (start) + ((end) - (start)) * t;
#define SMOOTH(start, end, t) LERP(start, end, (float)(-2.0f*pow((t), 3) + 3.0f*pow((t), 2)))
#define SMOOTH_2(start, end, t) LERP(start, end, (float)(6.0f*pow((t), 5) - 15.0f*pow((t), 4) + 10.0f*pow((t), 3)))

// Expects sample to be in [-1, 1]
inline VEC2 toDirection(float sample){
	sample += 1.0f;
	sample *= 0.5f;
	sample = clamp(sample);
	float theta = sample*2.0f*M_PI;
	return VEC2(cos(theta), sin(theta));
}

class PerlinNoise{
	VEC2 * gradient_;
	int width_, height_, max_grid_;
	float scale_;
public:
	PerlinNoise(const int& max_grid_in, const float& scale): scale_(scale){
		max_grid_ = max_grid_in;
		width_ = max_grid_ + 1;
		height_ = max_grid_ + 1;
		srand(12345);
		gradient_ = new VEC2[width_*height_];
		for(int i = 0; i < width_*height_; i++){
			int x = (int)rand()%8;
			VEC2 direction;
			switch(x){
				case 0:
					direction = VEC2(1.0f, 0.0f);
					break;
				case 1:
					direction = VEC2(1.0f, 1.0f);
					break;
				case 2:
					direction = VEC2(0.0f, 1.0f);
					break;
				case 3:
					direction = VEC2(-1.0f, 1.0f);
					break;
				case 4:
					direction = VEC2(-1.0f, 0.0f);
					break;
				case 5:
					direction = VEC2(-1.0f, -1.0f);
					break;
				case 6:
					direction = VEC2(0.0f, -1.0f);
					break;
				case 7:
					direction = VEC2(1.0f, -1.0f);
					break;
			}
			gradient_[i] = direction;
		}
	}

	~PerlinNoise(){
		delete[] gradient_;
	}

	float noise(const float& u, const float& v){
		float u_scaled = scale_*u;
		float v_scaled = scale_*v;

		VEC2 q = VEC2(u_scaled, v_scaled); // our sample position
		int min_u = (int)floor(u_scaled);
		int max_u = (min_u + 1);
		int min_v = (int)floor(v_scaled);
		int max_v = (min_v + 1);

		// looped u and v so that if we exceed the max grid we just repeat when finding gradient_ vectors
		int lMin_u = min_u % max_grid_;
		int lMax_u = max_u % max_grid_;
		int lMin_v = min_v % max_grid_;
		int lMax_v = max_v % max_grid_;


		int i_1 = min_v*width_ + min_u;
		int i_2 = min_v*width_ + max_u;
		int i_3 = max_v*width_ + max_u;
		int i_4 = max_v*width_ + min_u;

		// Arrays go counterclockwise starting with bottom left
		VEC2 pos[4] = {VEC2(min_u, min_v), VEC2(max_u, min_v), VEC2(max_u, max_v), VEC2(min_u, max_v)};
		VEC2 grad[4] = {gradient_[lMin_v*width_ + lMin_u], gradient_[lMin_v*width_ + lMax_u],
						gradient_[lMax_v*width_ + lMax_u], gradient_[lMax_v*width_ + lMin_u]}; // gradient_s vectors
		VEC2 dist[4];
		for(int i = 0; i < 4; i++) 
			dist[i] = pos[i] - q;
		float dot[4];
		for(int i = 0; i < 4; i++) 
			dot[i] = grad[i].dot(dist[i]);

		// Interpolate
		float mid_top = SMOOTH_2(dot[3], dot[2], (u_scaled-min_u)/1.0f);
		float mid_bottom = SMOOTH_2(dot[0], dot[1], (u_scaled-min_u)/1.0f);
		return SMOOTH_2(mid_bottom, mid_top, (v_scaled-min_v)/1.0f);
	}

	float turbulence(const float& u, const float& v){
		float N = (int)floor(log(max_grid_) / log(2));
		float sample = 0.0f;
		for(int i = 1; i <= N; i++){
			sample += noise(pow(2, i)*u, pow(2, i)*v) / pow(2, i);
		}
		return sample;
	}
};

class EvolvingNoise{
	PerlinNoise * pNoise_1;
	PerlinNoise * pNoise_2;
	float evolution_vel_;
public:
	EvolvingNoise(const int& max_grid_in, const float& scale, const float& evolution_vel): evolution_vel_(evolution_vel/scale){
		pNoise_1 = new PerlinNoise(max_grid_in, scale);
		pNoise_2 = new PerlinNoise(max_grid_in, scale);
	}

	float noise(const float& real_time, const float& u, const float& v){
		float evolution = evolution_vel_ * real_time;
		//std::cout << "evolution: " << evolution << "u: " << u << " v: " << v << std::endl;
		float sample = pNoise_1->noise(u + evolution, v);
		sample += pNoise_2->noise(u - evolution, v);
		sample /= 2;
		return sample;
	}
};

#endif // PerlinNoise_h