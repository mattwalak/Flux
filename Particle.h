#ifndef Particle_h
#define Particle_h

#include "SETTINGS.h"
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <limits>

#define IN_BOUNDS(x, xMax) (((x) >= 0) && ((x) < (xMax)))
#define COORD_IN_BOUNDS(x, y, xMax, yMax) (IN_BOUNDS((x), (xMax)) && IN_BOUNDS((y), (yMax)))
#define MAX_STEPS 1000000


/*

Notes to self:
write a method with a similar function to tracePath (Traces entire path of particle until off screen)
but put it in the tracer class instead (That's where the perlin noise object is for velocity)
ur prob gonna need to redesign this whole thing on paper because this is getting more complicated
that you original guessed
*/




/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class Particle{
	VEC2 pos_; // Real position/velocity
	VEC2 sim_pos_, sim_vel_; // Simulation position/velocity
	float start_time_, trail_length_, radius_, speed_;
	VEC3 color_;
public:
	void init(	const VEC2& pos, 
				const VEC3& color, 
				float start_time, 
				float speed,
				float trail_length, 
				float radius){
		pos_ = pos;
		color_ = color;
		start_time_ = start_time;
		speed_ = speed;
		trail_length_ = trail_length;
		radius_ = radius;
	}

	void setSimVel(const VEC2& vel_in) { sim_vel_ = vel_in; }
	VEC2 getSimPos() const { return sim_pos_; }
	float getStartTime() const { return start_time_; }
	VEC2 getPos() const { return pos_; }
	float getSpeed() const { return speed_; }
	float getTrailLength() const { return trail_length_; }
	VEC3 getColor() const { return color_; }
	float getRadius() const { return radius_; }

	// Advances the simulation
	// if we cross a pixel grid line durring this step, return true and set newPixel to coordinates of new pixel
	// else return false
	// newPixel is assumed to be a valid pointer to a VEC2
	bool step(const float& dt, VEC2 * newPixel){
		VEC2 nextPos = VEC2(sim_pos_[0] + sim_vel_[0]*dt, sim_pos_[1] + sim_vel_[1]*dt);
		if(floor(nextPos[0]) != floor(sim_pos_[0]) || floor(nextPos[1]) != floor(sim_pos_[1])){
			*newPixel = VEC2(floor(nextPos[0]), floor(nextPos[1]));
			sim_pos_ = nextPos;
			return true;
		}else{
			sim_pos_ = nextPos;
			return false;
		}
	}

	// Rewrite to account for stroke
	bool inBounds(int width, int height){
		return (sim_pos_[0] < width && sim_pos_[0] > 0.0f && sim_pos_[1] < height && sim_pos_[1] > 0.0f);
	}

	void resetSim(){
		sim_pos_ = pos_;
		sim_vel_ = VEC2(0.0f, 0.0f);
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class Tracer{
	// Constant properties
	int width_, height_;
	float granularity_;
	std::vector<Particle> particles_;
	EvolvingNoise * eNoise_;

	// Normalizes everything [0, 1]
	// assumes flux_out has already been allocated
	void normalizeFlux(VEC3 * toNormalize){
		float min = INFINITY;
		float max = -INFINITY;
		for(int i = 0; i < width_*height_; i++){
			float norm = toNormalize[i].norm();
			if(norm < min)
				min = norm;
			if(norm > max)
				max = norm;
		}
		float span = log(max + 1.0f) - log(min + 1.0f);
		for(int i = 0; i < width_*height_; i++){
			float norm = toNormalize[i].norm();
			toNormalize[i] = toNormalize[i].normalized() * clamp(log(norm + 1.0f)/span);
		}
	}

	// Assumes pixels have already been allocated
	// assumes norm_flux is already normalized
	void toImage(VEC3 * norm_flux, float * pixels){
		for(int y = 0; y < height_; y++){
			for(int x = 0; x < width_; x++){
				pixels[3*(y*width_ + x) + 0] = norm_flux[y*width_ + x][0] * 255.0f;
				pixels[3*(y*width_ + x) + 1] = norm_flux[y*width_ + x][1] * 255.0f;
				pixels[3*(y*width_ + x) + 2] = norm_flux[y*width_ + x][2] * 255.0f;
			}
		}
	}

	bool simParticle(float real_time, VEC3 * flux, Particle& p) const{
		p.resetSim();
		VEC2 pixel;
		float advance_time = real_time - p.getStartTime();
		if(advance_time < 0)
			return true;

		// find beginning and end of line
		int end_step = (int)(advance_time * p.getSpeed()/granularity_);
		int start_step = end_step - (int)(p.getTrailLength() * p.getSpeed()/granularity_);
		start_step = (start_step < 0 ? 0 : start_step);

		// Run steps
		for(int step = 0; step < end_step; step++){
			VEC2 pos = p.getSimPos();
			float sample = eNoise_->noise(real_time, pos[0]/width_, pos[1]/width_); // yes divide by width both times cause width will be larger
			p.setSimVel(toDirection(sample));
			bool result = p.step(granularity_, &pixel);

			if(step == start_step && !p.inBounds(width_, height_)){
				// We are already off the screen by the start of the tail -> this particle is done
				return false;
			}

			if(step >= start_step){
				// OPTION 1: Only update flux on pixel change
				/*if(result){
					// Draw code (Dot with radius stroke)
					int rad_bound = ceil(p.getRadius()) - 1;
					for(int dy = -rad_bound; dy <= rad_bound)


					int pix_x = (int)pixel[0];
					int pix_y = (int)pixel[1];
					if(COORD_IN_BOUNDS(pix_x, pix_y, width_, height_)){
						flux[pix_y*width_ + pix_x] += p.getColor();
					}
					
				}*/

				// OPTION 2: Always update flux
				float rad = p.getRadius();
				int rad_bound = ceil(rad);
				for(int dy = -rad_bound; dy <= rad_bound; dy++){
					for(int dx = -rad_bound; dx <= rad_bound; dx++){
						VEC2 real_pos = p.getSimPos();
						int pix_x = floor(real_pos[0])+dx;
						int pix_y = floor(real_pos[1])+dy;
						if(COORD_IN_BOUNDS(pix_x, pix_y, width_, height_)){
							VEC2 dipstick = VEC2(pix_x+0.5, pix_y+0.5);
							//float strength = 1.0f - (real_pos-dipstick).norm()/rad;
							//strength = (strength < 0 ? 0 : strength);
							flux[pix_y*width_ + pix_x] += p.getColor();
						}
					}
				}
			}
		}
		return true;
	}
public:
	Tracer(				int width, 
						int height, 
						float noise_scale, 
						float evolution_vel): 
						width_(width), height_(height) {
		eNoise_ = new EvolvingNoise(100, noise_scale, evolution_vel);
		granularity_ = 1.0f; // Default
	}

	void setGranularity(float gran) { granularity_ = gran; }

	void addParticle(	const VEC2& pos, 
						const VEC3& color, 
						const float& start_time, 
						const float& speed,
						const float& trail_length, 
						const float& radius){
		Particle p;
		p.init(pos, color, start_time, speed, trail_length, radius);
		particles_.push_back(p);
	}

	void render(float real_time, float * pixels){
		// Make the flux buffer
		VEC3 * flux = new VEC3[width_*height_](); // Gotta put this on the heap its so large
		for(int i = 0; i < width_*height_; i++) {
			flux[i] = VEC3(0.0f, 0.0f, 0.0f);
		}

		
		// Run simulation for each particle and update flux
		std::cout << "\t" << particles_.size() << "\n";
		for(int i = 0; i < particles_.size(); i++){
			if(!simParticle(real_time, flux, particles_[i])){
				particles_.erase(particles_.begin()+i);
				i--;
			}
		}
		/*
		// Just noise
		for(int y = 0; y < height_; y++){
			for(int x = 0; x < width_; x++){
				float u = (float)x/width_;
				float v = (float)y/width_;
				float sample = eNoise_->noise(real_time, u, v);
				sample /= 1.0f;
				sample += 0.5f;
				sample = clamp(sample);
				flux[y*width_ + x] = VEC3(1.0f, 1.0f, 1.0f) * sample;
			}
		}*/

		normalizeFlux(flux);
		toImage(flux, pixels);
		delete[] flux;
	}		
};

#endif // Particle_h