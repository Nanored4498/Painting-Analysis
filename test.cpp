#include "stb_image.h"
#include "lines.h"

int main() {
	int W, H, C;
	unsigned char* im = stbi_load("ims/SK-A-182.jpg", &W, &H, &C, 3);
	// unsigned char* im = stbi_load("ims/hooch1.jpg", &W, &H, &C, 3);
	// unsigned char* im = stbi_load("bilateral.png", &W, &H, &C, 3);
	PA::ProblemData* data = PA::applySobel(im, W, H);
	PA::save_sobel("sobel.png", data);
	return 1;
}