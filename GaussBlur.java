package iprocess;

public class GaussBlur {
	double[] kernel;

	public static void main(String[] args) throws Exception {
		GaussBlur shit = new GaussBlur();
		shit.setGaussianKernel(5, 1.6);
		ImageHandling h = new ImageHandling();
		double[][] data = h.getDataArray();
		double[][] newData = shit.convolve1D(data);
		h.saveArrayIntoImage(newData, "C://Users//STEVEN-PC//Desktop//img.jpg");
	}

	// sets kernel based on size of box (box must be odd); gaussian function
	public void setGaussianKernel(int box, double sigma) {
		kernel = new double[box];
		int radius = box / 2;
		int i = 0;
		for (int a = -radius; a <= radius; a++){
			kernel[i] = getGaussian(a, sigma);
			i++;
		}
	}

	// get values from gaussian function
	public double getGaussian(int index, double sigma) {
		return (1 / (Math.sqrt(2 * Math.PI) * (sigma))
				* Math.pow(Math.E, -Math.pow(index, 2) / (2 * Math.pow(sigma, 2))));
	}

	// passes 1D kernel in horizontal and vertical directions, due to separable
	// nature of gaussian
	public double[][] convolve1D(double[][] data) {
		int kernelSize = kernel.length / 2;
		int horiz = data.length;
		int vert = data[0].length;

		double[][] upd = new double[horiz][vert];
		// horizontal convolution 1D
		for (int i = 0; i < horiz; i++) {
			for (int j = 0; j < vert; j++) {
				double sum = 0;
				int index = 0;
				for (int k = -kernelSize; k <= kernelSize; k++) {
					int curr = i + k;
					if (curr < 0) {
						for (int m = 0; m < Math.abs(curr); m++) {
							sum += data[0][j] * kernel[index];
							index++;
						}
						k += Math.abs(curr);
					} else if (curr >= horiz) {
						sum += data[horiz - 1][j] * kernel[index];
						index++;
					} else { // 0 <= curr <= horiz
						sum += data[curr][j] * kernel[index];
						index++;
					}
				}

				upd[i][j] = sum;
				sum = 0;
				index = 0;
			}
		}

		double[][] aUpd = new double[horiz][vert];
		// vertical convolution 1D
		for (int i = 0; i < horiz; i++) {
			for (int j = 0; j < vert; j++) {
				double sum = 0;
				int index = 0;
				for (int k = -kernelSize; k <= kernelSize; k++) {
					int curr = j + k;
					if (curr < 0) {
						for (int m = 0; m < Math.abs(curr); m++) {
							sum += upd[i][0] * kernel[index];
							index++;
						}
						k += Math.abs(curr);
					} else if (curr >= vert) {
						sum += upd[i][vert - 1] * kernel[index];
						index++;
					} else { // 0 <= curr <= vert
						sum += upd[i][curr] * kernel[index];
						index++;
					}
				}

				aUpd[i][j] = sum;
				sum = 0;
				index = 0;
			}
		}

		return aUpd;
	}
	
	//extra functions i made 
	/**
	// Gaussian convolution using separable kernels
		public double[][] convGauss(double[][] intensities, double sigma) {
			int width = intensities.length;
			int height = intensities[0].length;

			double[] kernel = getKernel(sigma);

			double[][] convolvedImage = new double[width][height];

			convolvedImage = convolve1D(intensities, kernel);

			return convolvedImage;
		}

		// 1D Convolution in both directions
		public double[][] convolve1D(double[][] img, double[] kernel) {
			int width = img.length;
			int height = img[0].length;

			double[][] newimg = new double[width][height];

			int kernelRadius = kernel.length / 2;

			double[][] horzimg = new double[width][height];

			// convolve horizontal direction
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					horzimg[i][j] = convolvePixel(img, kernel, kernelRadius, width, height, i, j, true);
				}
			}

			// convolve vertical direction
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					newimg[i][j] = convolvePixel(horzimg, kernel, kernelRadius, width, height, i, j, false);
				}
			}

			return newimg;
		}

		// apply kernel at pixel (x,y) to obtain double val horizontal
		public double convolvePixel(double[][] img, double[] kernel, int radius, int width, int height, int x, int y,
				boolean horizontal) {
			double sum = 0;

			int count = 0;
			if (horizontal) {
				for (int i = x - radius; i < x + radius; i++) {

					// edge cases - wrap
					if (i < 0)
						sum += kernel[count] * img[width - 1 + i][y];
					else if (i >= width)
						sum += kernel[count] * img[width - i][y];
					else
						sum += kernel[count] * img[i][y];
				}
			} else { // vertical
				for (int j = y - radius; j < y + radius; j++) {

					// edge cases - wrap
					if (j < 0)
						sum += kernel[count] * img[x][height - 1 + j];
					else if (j >= width)
						sum += kernel[count] * img[x][height - j];
					else
						sum += kernel[count] * img[x][j];
				}
			}

			return sum;
		}

		// returns 1D Gaussian Kernel (since separable)
		public double[] getKernel(double sigma) {
			double[] arr = new double[2 * KERNEL_RADIUS + 1];
			int a = 0;
			for (int i = -KERNEL_RADIUS; i < KERNEL_RADIUS; i++) {
				double g = (double) (1 / (Math.sqrt(2 * Math.PI) * (sigma))
						* Math.pow(Math.E, -Math.pow(i, 2) / (2 * Math.pow(sigma, 2))));
				arr[a] = g;
				a++;
			}

			return arr;
		}
		**/

}
