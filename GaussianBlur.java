package iprocess;

import ij.ImagePlus;
import ij.io.OpenDialog;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;

public class GaussianBlur {

	/*
	 * Creates 1D Gaussian matrix to convolve image in x and y directions
	 */
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub

		String dir = findDir();
		ImagePlus im = new ImagePlus(dir);
		ImageProcessor imp = im.getProcessor();
		BufferedImage img = imp.getBufferedImage();

		for (int i = 0; i < 1; i++)
			img = gaussianBlur(img, 1, 2);

		File f = new File(dir);
		ImageIO.write(img, "png", f);
		System.out.println("SAVED");
	}

	public static BufferedImage gaussianBlur(BufferedImage img, double sigma, int box) {
		double[] kernel = getKernel(sigma, box);
		BufferedImage finalImage = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_RGB);

		// convolve each pixel

		// horizontal convolution
		for (int i = 0; i < img.getWidth() - 1; i++) {
			for (int j = 0; j < img.getHeight() - 1; j++) {
				int[] pixelrgb = conPixelHorizontal(img, kernel, box, i, j);
				int rgb = (pixelrgb[0] & 0x00FF) << 16 | (pixelrgb[1] & 0X00FF) << 8 | (pixelrgb[2] & 0X00FF);
				finalImage.setRGB(i, j, rgb);
			}
		}

		// vertical convolution
		BufferedImage finalImage2 = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < finalImage.getWidth() - 1; i++) {
			for (int j = 0; j < finalImage.getHeight() - 1; j++) {
				int[] pixelrgb = conPixelVertical(finalImage, kernel, box, i, j);
				int rgb = (pixelrgb[0] & 0x00FF) << 16 | (pixelrgb[1] & 0X00FF) << 8 | (pixelrgb[2] & 0X00FF);
				finalImage2.setRGB(i, j, rgb);
			}
		}

		return finalImage2;
	}

	public static int[] conPixelHorizontal(BufferedImage img, double[] kernel, int box, int x, int y) {
		int j = 0;

		double r = 0;
		double g = 0;
		double b = 0;
		for (int i = x - box; i <= x + box; i++) {
			double scaleVal = kernel[j];

			Color color;
			if (i < 0)
				color = new Color(img.getRGB(0, y));
			else if (i > img.getWidth() - 1)
				color = new Color(img.getRGB(img.getWidth() - 1, y));
			else
				color = new Color(img.getRGB(i, y));

			r += color.getRed() * scaleVal;
			g += color.getGreen() * scaleVal;
			b += color.getBlue() * scaleVal;

			j++;
		}

		int[] arr = new int[3];
		arr[0] = (int) r;
		arr[1] = (int) g;
		arr[2] = (int) b;
		return arr;
	}

	public static int[] conPixelVertical(BufferedImage img, double[] kernel, int box, int x, int y) {

		int j = 0;

		double r = 0;
		double g = 0;
		double b = 0;
		for (int i = y - box; i <= y + box; i++) {
			double scaleVal = kernel[j];

			Color color;
			if (i < 0)
				color = new Color(img.getRGB(x, 0));
			else if (i > img.getHeight() - 1)
				color = new Color(img.getRGB(x, img.getHeight() - 1));
			else
				color = new Color(img.getRGB(x, i));

			r += color.getRed() * scaleVal;
			g += color.getGreen() * scaleVal;
			b += color.getBlue() * scaleVal;

			j++;
		}
		int[] arr = new int[3];
		arr[0] = (int) r;
		arr[1] = (int) g;
		arr[2] = (int) b;
		return arr;
	}

	public static double[] getKernel(double sigma, int box) { // gets 1D kernel
																// according to
																// gaussian
																// function
		double[] arr = new double[2 * box + 1];
		int a = 0;
		for (int i = -box; i <= box; i++) {
			double g = (1 / Math.sqrt(2 * Math.PI * Math.pow(sigma, 2)))
					* Math.pow(Math.E, (-Math.pow(i, 2) / (2 * Math.pow(sigma, 2))));
			arr[a] = g;
			a++;
		}

		return arr;

	}

	public static String findDir() { // opens file opening dialog
		JFileChooser chooser = new JFileChooser();
		OpenDialog od = new OpenDialog("OPEN FILE", null);
		String dir = od.getDirectory();
		dir = dir.replace('\\', '/'); // Windows safe
		if (!dir.endsWith("/"))
			dir += "/";
		dir += od.getFileName();
		System.out.println(dir);

		return dir;
	}
}