# Flux
Some fun examples! (All rendered at 30fps, 1920x1080)

1000 particles with stroke = 2 with field evolution controlled by a sin wave (Evolution cycles seamlessly):
https://youtu.be/SLv5fj7NLWU

100 new particles every frame, stroke = 1, evolution velocity = 0.25:


How to make your own flux fields:
Place all particle generation code in the "ADD PARTICLES" section in the main() function of perlin.cpp
Frames are rendered one by one as .ppm images into the oven/ directory
Use ./convert (bash script) to convert .ppm sequence to video
