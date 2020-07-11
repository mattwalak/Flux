#ifndef Particle_h
#define Particle_h

#include "SETTINGS.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits>

#define IN_BOUNDS(x, xMax) (((x) >= 0) && ((x) < (xMax)))
#define COORD_IN_BOUNDS(x, y, xMax, yMax) (IN_BOUNDS((x), (xMax)) && IN_BOUNDS((y), (yMax)))

class Particle{
	VEC2 pos_, vel_;
public:
	void init(const VEC2& pos, const VEC2& vel){
		pos_ = pos;
		vel_ = vel;
	}

	void setPos(const VEC2& pos) {pos_ = pos; }
	VEC2 getPos() const { return pos_; }
	void setVel(const VEC2& vel) {vel_ = vel; }
	VEC2 getVel() const { return vel_; }

	// Advances the simulation
	// if we cross a pixel grid line durring this step, return true and set newPixel to coordinates of new pixel
	// else return false
	// newPixel is assumed to be a valid pointer to a VEC2
	bool step(const float& dt, VEC2 * newPixel){
		VEC2 nextPos = VEC2(pos_[0] + vel_[0]*dt, pos_[1] + vel_[1]*dt);
		if(floor(nextPos[0]) != floor(pos_[0]) || floor(nextPos[1]) != floor(pos_[1])){
			*newPixel = VEC2(floor(nextPos[0]), floor(nextPos[1]));
			pos_ = nextPos;
			return true;
		}else{
			pos_ = nextPos;
			return false;
		}
	}
};

class Tracer{
	int width_, height_;
	float * flux_;
	float scale_;
	PerlinNoise * pNoise_;
	std::vector<Particle> particles_;

	float field_strength = 1.0f;
	float step_dt = 2.0f;
	int particle_radius = 2;
	float max_edge_dist = sqrt(2.0f*pow(particle_radius, 2));
public:
	Tracer(int width, int height, int numParticles, float scale): width_(width), height_(height), scale_(scale) {
		flux_ = new float[width_*height_]();
		pNoise_ = new PerlinNoise(10);

		/*
		for(int x = 0; x < width_; x++){
			for(int y = 0; y < height_; y++){
				for(int i = 0; i < 4; i++){
					for(int j = 0; j < 4; j++){
						Particle p;
						p.init(VEC2(x + i/4.0f, y + j/4.0f), VEC2(0.0f, 0.0f));
						particles_.push_back(p);
					}
				}
			}
		}*/

		for(int i = 0; i < numParticles; i++){
			Particle p;
			p.init(VEC2(rand()%width_, rand()%height), VEC2(0.0f, 0.0f));
			particles_.push_back(p);
		}
	}

	// Updates each particle with new velocity based on new position
	void updateVelocity(){
		for(int i = 0; i < particles_.size(); i++){
			VEC2 particle_pos = particles_[i].getPos();
			float u = scale_*particle_pos[0]/width_;
			float v = scale_*particle_pos[1]/height_;
			float sample = pNoise_->noise(u, v);
			VEC2 newVel = toDirection(sample);
			particles_[i].setVel(newVel);
		}
	}

	// Steps each particle forward and removes particle if off screen
	void step(){
		for(int i = 0; i < particles_.size(); i++){
			VEC2 pixel;
			bool result = particles_[i].step(step_dt, &pixel);
			VEC2 newPos = particles_[i].getPos();
			if(newPos[0] > width_ || newPos[0] < 0.0f || newPos[1] > height_ || newPos[1] < 0.0f){
				// Particle is not on screen
				particles_.erase(particles_.begin()+i);
				i--;
			}else{
				// Particle is still on screen (particle is 3x3 box lol)
				if(result){
					int pix_x = (int)pixel[0];
					int pix_y = (int)pixel[1];
					for(int x = -particle_radius; x <= particle_radius; x++){
						for(int y = -particle_radius; y <= particle_radius; y++){
							int coord_x = pix_x + x;
							int coord_y = pix_y + y;
							if(COORD_IN_BOUNDS(coord_x, coord_y, width_, height_)){
								float dist = sqrt(pow(abs(x), 2) + pow(abs(y), 2));
								float strength = 1.0f - dist/max_edge_dist;
								flux_[coord_y*width_ + coord_x] += strength;
							}
						}
					}

					// flux_[pix_y*width_ + pix_x] += 1.0f; // pixel increases flux!
				}
			}
		}
	}

	int particleCount(){
		return particles_.size();
	}

	void printStats(){
		std::cout << "Num particles = " << particles_.size() << std::endl;
		for(int i = 0; i < particles_.size(); i++){
			VEC2 pos = particles_[i].getPos();
			VEC2 vel = particles_[i].getVel();
			std::cout << i << "): pos = <" << pos[0] << ", " << pos[1] << "> vel = <" << vel[0] << ", " << vel[1] << ">\n";
		}
	}

	// Normalizes everything [0, 1]
	// assumes flux_out has already been allocated
	void getNormalizedFlux(float * flux_out){
		float min = INFINITY;
		float max = -INFINITY;
		for(int i = 0; i < width_*height_; i++){
			if(flux_[i] < min)
				min = flux_[i];
			if(flux_[i] > max)
				max = flux_[i];
		}

		float span = log(max + 1.0f) - log(min + 1.0f);
		for(int i = 0; i < width_*height_; i++){
			flux_out[i] = clamp(log(flux_[i] + 1.0f)/span);
		}
	}

	// Assumes pixels have already been allocated
	void toImage(float * pixels){
		float norm_flux[width_*height_];
		getNormalizedFlux(norm_flux);
		for(int y = 0; y < height_; y++){
			for(int x = 0; x < width_; x++){
				pixels[3*(y*width_ + x) + 0] = norm_flux[y*width_ + x] * 255.0f;
				pixels[3*(y*width_ + x) + 1] = norm_flux[y*width_ + x] * 255.0f;
				pixels[3*(y*width_ + x) + 2] = norm_flux[y*width_ + x] * 255.0f;
			}
		}
	}	
};

#endif // Particle_h