package iprocess;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;

import ij.ImagePlus;
import ij.io.OpenDialog;

public class SIFT {

	static int KERNEL_RADIUS = 2;
	static double BASE_SIGMA = Math.sqrt(2);
	static int NUM_SCALE = 5;
	static int NUM_OCTAVES = 4;
	static double START_SIGMA = 1.6;
 
	public static void main(String[] args) {
		// TODO Auto-generated method stub
	}

	// ** Difference of Gaussians **

	public ArrayList<ArrayList<double[][]>> createScaleSpace(double[][] img) {
		ArrayList<ArrayList<double[][]>> scaleSpace = new ArrayList<ArrayList<double[][]>>();

		double[][] currimg = img.clone();
		for (int i = 0; i < NUM_OCTAVES; i++) {
			ArrayList<double[][]> data = buildDOGPyramid(currimg);
			scaleSpace.add(data);
			currimg = subSample(currimg);
		}

		return scaleSpace;
	}

	// Builds DOG pyramid at certain octave
	public ArrayList<double[][]> buildDOGPyramid(double[][] img) {
		ArrayList dogPyrm = new ArrayList<double[][]>();
		ArrayList gaussianPyrm = new ArrayList<double[][]>();

		int width = img.length;
		int height = img[0].length;

		for (int i = 0; i < NUM_SCALE; i++) {
			double[][] newImg = convGauss(img, BASE_SIGMA);
			gaussianPyrm.add(newImg);
			img = newImg.clone();
		}

		return laplacianDifference(gaussianPyrm, width, height);

	}

	// Takes difference between L(x,y,o) - L(x,y,ko)
	public ArrayList<double[][]> laplacianDifference(ArrayList<double[][]> pyrm, int width, int height) {
		int size = pyrm.size();
		ArrayList<double[][]> dogPyrm = new ArrayList<double[][]>();
		for (int i = 0; i < size; i++){
			dogPyrm.add(new double[width][height]);
		}
		
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				for (int k = 0; k < size; k++) {
//					dogPyrm.
				}
			}
		}
		
		return null;
	}

	
	// halves image dimensions
	public static double[][] subSample(double[][] img) {
		int width = img.length;
		int height = img[0].length;

		double[][] halvedImg = new double[width / 2][height / 2];
		for (int i = 0; i < width; i += 2) {
			for (int j = 0; j < height; j += 2) {
				halvedImg[i / 2][j / 2] = img[i][j];
			}
		}

		return halvedImg;
	}

}
