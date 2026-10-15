// This will apply the sobel filter and return the PSNR between the golden sobel and the produced sobel
// sobelized image
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

#define SIZE	4096
#define INPUT_FILE	"input.grey"
#define OUTPUT_FILE	"output_sobel.grey"
#define GOLDEN_FILE	"golden.grey"

#define pow2(x) ({		\
	x * x;				\
})

// This is a fast approximation to log2()
// Y = C[0]*F*F*F + C[1]*F*F + C[2]*F + C[3] + E;
float log2f_approx(float X) {
  float Y, F;
  int E;
  F = frexpf(fabsf(X), &E);
  Y = 1.23149591368684f;
  Y *= F;
  Y += -4.11852516267426f;
  Y *= F;
  Y += 6.02197014179219f;
  Y *= F;
  Y += -3.13396450166353f;
  Y += E;
  return(Y);
}

double sobel(unsigned char *input, unsigned char *output, unsigned char *golden);

/* The arrays holding the input image, the output image and the output used *
 * as golden standard. The luminosity (intensity) of each pixel in the      *
 * grayscale image is represented by a value between 0 and 255 (an unsigned *
 * character). The arrays (and the files) contain these values in row-major *
 * order (element after element within each row and row after row. 			*/
unsigned char  input[SIZE*SIZE], output[SIZE*SIZE];
/*const*/ unsigned char golden[SIZE*SIZE];

/* Implement a 2D convolution of the matrix with the operator */
/* posy and posx correspond to the vertical and horizontal disposition of the *
 * pixel we process in the original image, input is the input image and       *
 * operator the operator we apply (horizontal or vertical). The function ret. *
 * value is the convolution of the operator with the neighboring pixels of the*
 * pixel we process.														  */

/* The main computational function of the program. The input, output and *
 * golden arguments are pointers to the arrays used to store the input   *
 * image, the output produced by the algorithm and the output used as    *
 * golden standard for the comparisons.									 */
double sobel(unsigned char *input, unsigned char *output, unsigned char *golden)
{
	register long double PSNR = 0;
	int i, j;
	unsigned int p;
	int res;
	struct timespec  tv1, tv2;
	FILE *f_in, *f_out, *f_golden;

	/* The first and last row of the output array, as well as the first  *
     * and last element of each column are not going to be filled by the *
     * algorithm, therefore make sure to initialize them with 0s.		 */
	memset(output, 0, SIZE*sizeof(unsigned char));
	memset(&output[SIZE*(SIZE-1)], 0, SIZE*sizeof(unsigned char));

	/* Open the input, output, golden files, read the input and golden    *
     * and store them to the corresponding arrays.						  */
	f_in = fopen(INPUT_FILE, "r");
	if (f_in == NULL) {
		printf("File " INPUT_FILE " not found\n");
		exit(1);
	}

	f_out = fopen(OUTPUT_FILE, "wb");
	if (f_out == NULL) {
		printf("File " OUTPUT_FILE " could not be created\n");
		fclose(f_in);
		exit(1);
	}  

	f_golden = fopen(GOLDEN_FILE, "r");
	if (f_golden == NULL) {
		printf("File " GOLDEN_FILE " not found\n");
		fclose(f_in);
		fclose(f_out);
		exit(1);
	}    

	fread(input, sizeof(unsigned char), SIZE*SIZE, f_in);
	fread(golden, sizeof(unsigned char), SIZE*SIZE, f_golden);
	fclose(f_in);
	fclose(f_golden);

	/* This is the main computation. Get the starting time. */
	clock_gettime(CLOCK_MONOTONIC_RAW, &tv1);
	/* For each pixel of the output image */

	int iSize = SIZE;
	for (i=1; i<SIZE-1; i+=1) {

		output[iSize] = 0;
		output[iSize + SIZE - 1] = 0;
		for (j=1; j<SIZE-1; j+=8) {
			/* Apply the sobel filter and calculate the magnitude *
			 * of the derivative.	*/	

			unsigned char sobel1 = input[iSize - SIZE + j-1];
			unsigned char sobel2 = input[iSize - SIZE + j  ];
			unsigned char sobel3 = input[iSize - SIZE + j+1];
			unsigned char sobel4 = input[ iSize + j-1];
			unsigned char sobel5 = input[ iSize + j  ];
			unsigned char sobel6 = input[ iSize + j+1];
			unsigned char sobel7 = input[iSize + SIZE + j-1];
			unsigned char sobel8 = input[iSize + SIZE + j  ];
			unsigned char sobel9 = input[iSize + SIZE + j+1];
			unsigned char sobel10 = input[iSize - SIZE + (j+2)];
			unsigned char sobel11 = input[iSize + (j+2)];
			unsigned char sobel12 = input[iSize + SIZE + (j+2)];
			unsigned char sobel13 = input[iSize - SIZE + (j+3)];
			unsigned char sobel14 = input[iSize + (j+3)];
			unsigned char sobel15 = input[iSize + SIZE + (j+3)];
			unsigned char sobel16 = input[iSize - SIZE + (j+4)];
			unsigned char sobel17 = input[iSize + (j+4)];
			unsigned char sobel18 = input[iSize + SIZE + (j+4)];

			unsigned char sobel19 = input[iSize - SIZE + (j+5)];
			unsigned char sobel20 = input[iSize + (j+5)];
			unsigned char sobel21 = input[iSize + SIZE + (j+5)];
			unsigned char sobel22 = input[iSize - SIZE + (j+6)];
			unsigned char sobel23 = input[iSize + (j+6)];
			unsigned char sobel24 = input[iSize + SIZE + (j+6)];
			unsigned char sobel25 = input[iSize - SIZE + (j+7)];
			unsigned char sobel26 = input[iSize + (j+7)];
			unsigned char sobel27 = input[iSize + SIZE + (j+7)];
			unsigned char sobel28 = input[iSize - SIZE + (j+8)];
			unsigned char sobel29 = input[iSize + (j+8)];
			unsigned char sobel30 = input[iSize + SIZE + (j+8)];

    		int conv1H = 
			- sobel1 
			+ sobel3 
			- sobel4 - sobel4  
			+ sobel6 + sobel6  
			- sobel7  
			+ sobel9 ; 	
			
			int conv1V = 
			+ sobel1  
			+ sobel2 + sobel2 
			+ sobel3 
			- sobel7  
			- sobel8 - sobel8 
			- sobel9 ;

			p = pow2(conv1H) + pow2(conv1V);
			res = (int)sqrt(p);

			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + j] = (unsigned char)res;
			else
				output[iSize + j] = 255;

			// t = pow2((output[iSize+j] - golden[iSize+j]));
			PSNR += pow2((output[iSize+(j)] - golden[iSize+(j)]));

			int conv2H = 
			- sobel2 
			+ sobel10 
			- sobel5 - sobel5 
			+ sobel11 + sobel11 
			- sobel8  
			+ sobel12;

			int conv2V = 
			+ sobel2  
			+ sobel3 + sobel3 
			+ sobel10 
			- sobel8 
			- sobel9 - sobel9 
			- sobel12;

			p = pow2(conv2H) + pow2(conv2V);
			res = (int)sqrt(p);
			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+1)] = (unsigned char)res;
			else
				output[iSize + (j+1)] = 255;

			// t = pow2((output[iSize+(j+1)] - golden[iSize+(j+1)]));
			PSNR += pow2((output[iSize+(j+1)] - golden[iSize+(j+1)]));

			int conv3H = 
			- sobel3 
			+ sobel13 
			- sobel6 - sobel6 
			+ sobel14 + sobel14 
			- sobel9  
			+ sobel15;

			int conv3V = 
			+ sobel3  
			+ sobel10 + sobel10 
			+ sobel13 
			- sobel9 
			- sobel12 - sobel12 
			- sobel15;

			p = pow2(conv3H) + pow2(conv3V);
			res = (int)sqrt(p);
			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+2)] = (unsigned char)res;
			else
				output[iSize + (j+2)] = 255;
			

			// t = pow2((output[iSize+(j+2)] - golden[iSize+(j+2)]));
			PSNR += pow2((output[iSize+(j+2)] - golden[iSize+(j+2)]));

			int conv4H = 
			- sobel10 
			+ sobel16 
			- sobel11 - sobel11 
			+ sobel17 + sobel17 
			- sobel12  
			+ sobel18;

			int conv4V = 
			+ sobel10  
			+ sobel13 + sobel13 
			+ sobel16 
			- sobel12 
			- sobel15 - sobel15 
			- sobel18;

			p = pow2(conv4H) + pow2(conv4V);
			res = (int)sqrt(p);
			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+3)] = (unsigned char)res;
			else
				output[iSize + (j+3)] = 255;

			// t = pow2((output[iSize+(j+3)] - golden[iSize+(j+3)]));
			PSNR += pow2((output[iSize+(j+3)] - golden[iSize+(j+3)]));

			int conv5H = 
			- sobel13 
			+ sobel19 
			- sobel14 - sobel14  
			+ sobel20 + sobel20  
			- sobel15  
			+ sobel21 ; 	

			int conv5V = 
			+ sobel13  
			+ sobel16 + sobel16 
			+ sobel19 
			- sobel18  
			- sobel21 - sobel21 
			- sobel24 ;

			p = pow2(conv5H) + pow2(conv5V);
			res = (int)sqrt(p);

			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+4)] = (unsigned char)res;
			else
				output[iSize + (j+4)] = 255;

			// t = pow2((output[iSize+j] - golden[iSize+j]));
			PSNR += pow2((output[iSize+(j+4)] - golden[iSize+(j+4)]));

			int conv6H = 
			- sobel16 
			+ sobel22 
			- sobel17 - sobel17  
			+ sobel23 + sobel23  
			- sobel18  
			+ sobel24 ; 	
			
			int conv6V = 
			+ sobel16  
			+ sobel19 + sobel19 
			+ sobel22 
			- sobel18  
			- sobel21 - sobel21 
			- sobel24 ;

			p = pow2(conv6H) + pow2(conv6V);
			res = (int)sqrt(p);

			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+5)] = (unsigned char)res;
			else
				output[iSize + (j+5)] = 255;

			// t = pow2((output[iSize+j] - golden[iSize+j]));
			PSNR += pow2((output[iSize+(j+5)] - golden[iSize+(j+5)]));

			int conv7H = 
			- sobel19 
			+ sobel25 
			- sobel20 - sobel20  
			+ sobel26 + sobel26  
			- sobel21  
			+ sobel27 ; 	
			
			int conv7V = 
			+ sobel19  
			+ sobel22 + sobel22 
			+ sobel25 
			- sobel21  
			- sobel24 - sobel24 
			- sobel27 ;

			p = pow2(conv7H) + pow2(conv7V);
			res = (int)sqrt(p);

			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+6)] = (unsigned char)res;
			else
				output[iSize + (j+6)] = 255;

			// t = pow2((output[iSize+j] - golden[iSize+j]));
			PSNR += pow2((output[iSize+(j+6)] - golden[iSize+(j+6)]));

			int conv8H = 
			- sobel22  
			+ sobel28 
			- sobel23 - sobel23
			+ sobel29 + sobel29 
			- sobel24
			+ sobel30; 

			int conv8V = 
			+ sobel22  
			+ sobel25 + sobel25 
			+ sobel28 
			- sobel24  
			- sobel27 - sobel27 
			- sobel30;

			p = pow2(conv8H) + pow2(conv8V);
			res = (int)sqrt(p);

			/* If the resulting value is greater than 255, clip it *
			 * to 255.											   */
			if(res < 255)
				output[iSize + (j+7)] = (unsigned char)res;
			else
				output[iSize + (j+7)] = 255;

			// t = pow2((output[iSize+j] - golden[iSize+j]));
			PSNR += pow2((output[iSize+(j+7)] - golden[iSize+(j+7)]));
		}
		iSize += SIZE;
	}

	/* Now run through the output and the golden output to calculate *
	 * the MSE and then the PSNR.									 */

	PSNR /= (double)(SIZE*SIZE);
	const double tmp2 = log2f(65536/PSNR);
	PSNR =  + tmp2 
			+ tmp2
			+ tmp2 
			+ tmp2
			+ tmp2
			+ tmp2
			+ tmp2 
			+ tmp2
			+ tmp2
			+ tmp2;

	/* This is the end of the main computation. Take the end time,  *
	 * calculate the duration of the computation and report it. 	*/
	clock_gettime(CLOCK_MONOTONIC_RAW, &tv2);

	printf ("Total time = %10g seconds\n",
			(double) (tv2.tv_nsec - tv1.tv_nsec) / 1000000000.0 +
			(double) (tv2.tv_sec - tv1.tv_sec));


	/* Write the output file */
	fwrite(output, sizeof(unsigned char), SIZE*SIZE, f_out);
	fclose(f_out);
  
	return PSNR;
}


int main(int argc, char* argv[])
{
	double PSNR;
	PSNR = sobel(input, output, golden);
	printf("PSNR of original Sobel and computed Sobel image: %g\n", PSNR);
	printf("A visualization of the sobel filter can be found at " OUTPUT_FILE ", or you can run 'make image' to get the jpg\n");

	return 0;
}

