#ifndef COLOR_h
#define COLOR_h

#include<cmath> // Needed for fmod()
#include <vector>

/*
 * H(Hue): 0 - 360 degree (integer)
 * S(Saturation): 0 - 1.00 (double)
 * V(Value): 0 - 1.00 (double)
 * 
 * output[3]: Output, array size 3, int
 */
// Credit: https://gist.github.com/kuathadianto/200148f53616cbd226d993b400214a7f
void HSVtoRGB(int H, double S, double V, float output[3]) {
	double C = S * V;
	double X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	double m = V - C;
	double Rs, Gs, Bs;

	if(H >= 0 && H < 60) {
		Rs = C;
		Gs = X;
		Bs = 0;	
	}
	else if(H >= 60 && H < 120) {	
		Rs = X;
		Gs = C;
		Bs = 0;	
	}
	else if(H >= 120 && H < 180) {
		Rs = 0;
		Gs = C;
		Bs = X;	
	}
	else if(H >= 180 && H < 240) {
		Rs = 0;
		Gs = X;
		Bs = C;	
	}
	else if(H >= 240 && H < 300) {
		Rs = X;
		Gs = 0;
		Bs = C;	
	}
	else {
		Rs = C;
		Gs = 0;
		Bs = X;	
	}
	
	output[0] = (Rs + m);
	output[1] = (Gs + m);
	output[2] = (Bs + m);
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

#endif // COLOR_h
