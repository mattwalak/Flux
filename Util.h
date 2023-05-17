#ifndef Util_h
#define Util_h

#include <iostream>
#include <string>

float clamp(const float& t){
  if(t > 1.0f)
    return 1.0f;
  else if(t < 0.0f)
    return 0.0f;
  else
    return t;
}

void writePPM(const std::string& filename, int xRes, int yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    std::cout << " Could not open file \"" << filename.c_str() << "\" for writing." << std::endl;
    std::cout << " Make sure you're not trying to write from a weird location or with a " << std::endl;
    std::cout << " strange filename. Bailing ... " << std::endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}

void readPPM(const std::string& filename, int xRes, int yRes, float*& values)
{
  // try to open the file
  FILE *fp;
  fp = fopen(filename.c_str(), "rb");
  if (fp == NULL)
  {
    std::cout << " Could not open file \"" << filename.c_str() << "\" for reading." << std::endl;
    std::cout << " Make sure you're not trying to read from a weird location or with a " << std::endl;
    std::cout << " strange filename. Bailing ... " << std::endl;
    exit(0);
  }

  // get the dimensions
  fscanf(fp, "P6\n%d %d\n255\n", &xRes, &yRes);
  int totalCells = xRes * yRes;

  // grab the pixel values
  unsigned char* pixels = new unsigned char[3 * totalCells];
  fread(pixels, 1, totalCells * 3, fp);

  // copy to a nicer data type
  values = new float[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    values[i] = pixels[i];

  // clean up
  delete[] pixels;
  fclose(fp);
  std::cout << " Read in file " << filename.c_str() << std::endl;
}


#endif // Util_h