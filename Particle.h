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
	VEC2 pos_, vel_;
	VEC3 color_;
	float trail_length_sec_, steps_per_sec_, start_time_;

	int particle_radius_ = 2; // Construct this somewhere please
public:
	void init(const VEC2& pos, const VEC2& vel, const VEC3& color, float trail_length_sec, float steps_per_sec, float start_time){
		pos_ = pos;
		vel_ = vel;
		color_ = color;
		trail_length_sec_ = trail_length_sec;
		steps_per_sec_ = steps_per_sec;
		start_time_ = start_time;
	}

	void setPos(const VEC2& pos) {pos_ = pos; }
	VEC2 getPos() const { return pos_; }
	void setVel(const VEC2& vel) {vel_ = vel; }
	VEC2 getVel() const { return vel_; }
	void setColor(const VEC3& color) {color_ = color; }
	VEC3 getColor() const { return color_; }

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

	// Traces full flux path and adds it to the buffer
	// Returns true to keep on screen for next frame
	// Return false to destroy before next frame
	bool tracePath(const float& real_time, int width, int height, const float& dt, VEC3 * flux){
		if(real_time < start_time_)
			return true; // Perhaps you're not born yet

		// advance steps to account for real_time offset
		VEC2 newPixel;
		int step_advance = (int)(real_time - start_time_)*steps_per_sec_;
		for(int i = 0; i < step_advance; i++){
			step(dt, &newPixel);
		}

		if(!inBounds(width, height))
			return false; // Destroy next frame! (Nothing left to draw)

		while(inBounds(width, height)){
			if(step(dt, &newPixel)){
				// Circle drawing code
				int pix_x = (int)newPixel[0];
				int pix_y = (int)newPixel[1];
				for(int x = -particle_radius_; x <= particle_radius_; x++){
					for(int y = -particle_radius_; y <= particle_radius_; y++){
						int coord_x = pix_x + x;
						int coord_y = pix_y + y;
						if(COORD_IN_BOUNDS(coord_x, coord_y, width, height)){
							float dist = sqrt(pow(abs(x), 2) + pow(abs(y), 2));
							if(dist <= particle_radius_){
								float strength = 1.0f - dist/particle_radius_;
								flux[coord_y*width + coord_x] += strength * color_;
							}
						}
					}
				}
			}
		}
		return true;
	}

	// Rewrite to account for stroke
	bool inBounds(int width, int height){
		return (pos_[0] < width && pos_[0] > 0.0f && pos_[1] < height && pos_[1] > 0.0f);
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
class SimParticle{
	float trail_length_sec_;
	float velocity_; // Steps per second
	float start_time_; // When this particle was born
	int step_count_;

	Particle initial_; // Initial position of this particle (does not change)
	Particle sim_; // Simulation position (Changes each step)
public:
	// All getters and setters all refer to the sim_ particle
	void setPos(const VEC2& pos) {sim_.setPos(pos); }
	VEC2 getPos() const { return sim_.getPos(); }
	void setVel(const VEC2& vel) {sim_.setVel(vel); }
	VEC2 getVel() const { return sim_.getVel(); }
	void setColor(const VEC3& color) {sim_.setColor(color); }
	VEC3 getColor() const { return sim_.getColor(); }

	// Init() is the only function that changes initial_
	void init(const VEC2& pos, const VEC2& vel, const VEC3& color, const float& trail_length_sec, const float& velocity, const float& current_time){
		initial_.setPos(pos);
		initial_.setVel(vel);
		initial_.setColor(color);
		sim_.setPos(pos);
		sim_.setVel(vel);
		sim_.setColor(color);
		trail_length_sec_ = trail_length_sec;
		velocity_ = velocity;
		step_count_ = 0;
		start_time_ = current_time;
	}

	// Reset for next frame
	void reset(){
		sim_ = initial_;
		step_count_ = 0;
	}

	// Advances the simulation
	// if we cross a pixel grid line durring this step, return true and set newPixel to coordinates of new pixel
	// else return false
	// newPixel is assumed to be a valid pointer to a VEC2
	// Refers to sim_... obviously...
	bool step(float real_t, float dt, VEC2 * newPixel){
		return sim_.step(dt, newPixel);
	}

};*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class Tracer{
	// Constant properties
	int width_, height_;
	EvolvingNoise * eNoise_;

	// Stores all active particles
	std::vector<Particle> particles_;

	// Stores only the on-screen particles in a sim group
	std::vector<Particle> sim_group_; 

	// Flux buffer
	VEC3 * flux_;

	// Animatable properties
	float real_time_;
	int steps_per_sec_ = 100;
	float trail_length_sec_ = 1.0f;
	float scale_;
	std::vector<VEC3> colors_;
	float step_dt_ = 2.0f;
	int particle_radius_ = 2;
	float max_edge_dist_ = sqrt(2.0f*pow(particle_radius_, 2));
public:
	// Setters

	Tracer(int width, int height, float scale, const std::vector<VEC3>& colors): width_(width), height_(height), scale_(scale) {
		flux_ = new VEC3[width_*height_]();
		eNoise_ = new EvolvingNoise(10);
		colors_ = colors;
	}

	// Note we are adding at real_time, not real_time_ (It is possible to add particles in the past / future)
	void addParticle(float real_time, float velocity, float trail_length_sec){
		int color_i = rand()%colors_.size();
		float pos_x = rand()%width_;
		float pos_y = rand()%height_;
		std::cout << "new particle at x: " << pos_x << ", y: " << pos_y << std::endl;
		VEC2 pos = VEC2(pos_x, pos_y);
		VEC2 vel = VEC2(0.0f, 0.0f);
		VEC3 color = colors_[color_i];
		Particle p;
		p.init(pos, vel, color, trail_length_sec, velocity, real_time);
		particles_.push_back(p);
	}

	// Updates each particle with new velocity based on new position
	void updateVelocity(float evolution){
		for(int i = 0; i < sim_group_.size(); i++){
			VEC2 particle_pos = sim_group_[i].getPos();
			float u = scale_*particle_pos[0]/width_;
			float v = scale_*particle_pos[1]/height_;
			//float sample = pNoise_->noise(u, v);
			float sample = eNoise_->noise(u, v, evolution);
			VEC2 newVel = toDirection(sample);
			sim_group_[i].setVel(newVel);
		}
	}

	// Steps each particle forward and removes particle if off screen
	void step(){
		for(int i = 0; i < sim_group_.size(); i++){
			Particle * p = &sim_group_[i];

			// Perform step
			VEC2 pixel;
			bool result = p->step(step_dt_, &pixel);

			// Check if in bounds
			if(!p->inBounds(width_, height_)){
				// Particle is not on screen
				sim_group_.erase(sim_group_.begin()+i);
				i--;
			}else{
				// Particle is still on screen
				if(result){
					// Circle drawing code
					int pix_x = (int)pixel[0];
					int pix_y = (int)pixel[1];
					for(int x = -particle_radius_; x <= particle_radius_; x++){
						for(int y = -particle_radius_; y <= particle_radius_; y++){
							int coord_x = pix_x + x;
							int coord_y = pix_y + y;
							if(COORD_IN_BOUNDS(coord_x, coord_y, width_, height_)){
								float dist = sqrt(pow(abs(x), 2) + pow(abs(y), 2));
								if(dist <= particle_radius_){
									float strength = 1.0f - dist/particle_radius_;
									flux_[coord_y*width_ + coord_x] += strength * p->getColor();
								}
							}
						}
					}
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
	void getNormalizedFlux(VEC3 * flux_out){
		float min = INFINITY;
		float max = -INFINITY;
		for(int i = 0; i < width_*height_; i++){
			float norm = flux_[i].norm();
			if(norm < min)
				min = norm;
			if(norm > max)
				max = norm;
		}
		float span = log(max + 1.0f) - log(min + 1.0f);
		for(int i = 0; i < width_*height_; i++){
			float norm = flux_[i].norm();
			flux_out[i] = flux_[i].normalized() * clamp(log(norm + 1.0f)/span);
		}
	}

	// Called every time before simulateFrame / toImage sequence
	void resetBuffers(){
		// Clear buffers
		for(int i = 0; i < width_*height_; i++) flux_[i] = VEC3(0.0f, 0.0f, 0.0f);
	}

	void setRealTime(float real_time){
		real_time_ = real_time;
	}

	// Run the simulation!
	void simulateFrame(float evolution){
		sim_group_ = particles_; // copy booiiii
		for(int i = 0; i < sim_group_.size(); i++){
			if(!sim_group_[i].tracePath(real_time_, width_, height_, step_dt_, flux_)){
				sim_group_.erase(sim_group_.begin()+i);
				i--;
			}
		}



		for(int step_num = 0; step_num < MAX_STEPS; step_num++){
			if(particleCount() == 0)
				return;
			if(step_num == 0)
				std::cout << "Running sim" << std::endl;

			updateVelocity(evolution);
			step();
		}
	}

	// Assumes pixels have already been allocated
	void toImage(float * pixels){
		//VEC3 norm_flux[1000*1000];
		VEC3 * norm_flux = new VEC3[width_*height_](); // I am so confused why I need this
		getNormalizedFlux(norm_flux);
		for(int y = 0; y < height_; y++){
			for(int x = 0; x < width_; x++){
				pixels[3*(y*width_ + x) + 0] = norm_flux[y*width_ + x][0] * 255.0f;
				pixels[3*(y*width_ + x) + 1] = norm_flux[y*width_ + x][1] * 255.0f;
				pixels[3*(y*width_ + x) + 2] = norm_flux[y*width_ + x][2] * 255.0f;
			}
		}
		delete[] norm_flux; // Like wtf
	}	
};

#endif // Particle_h