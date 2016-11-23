package iprocess;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.imageio.ImageIO;

import org.apache.commons.lang3.ArrayUtils;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import ij.ImagePlus;

public class EigenFaces {
	final static int DIMENSION = 320 * 243;
	
	public static void main(String[] args) throws Exception {
		File[] a = getFiles("C://Users//STEVEN-PC//Desktop//centerimgs");
		Matrix dataMatrix = getDataMatrix(a);
		Matrix meanVector = getMeanVector(dataMatrix);
		Matrix differenceMatrix = getDifferenceMatrix(meanVector, dataMatrix); //A Matrix
		Matrix covEigVec = getEigenface(differenceMatrix);
		Matrix omega = getWeights(covEigVec, differenceMatrix);
//		Matrix test = getMatrixFromFile("C://Users//STEVEN-PC//Desktop//yalefaces//subject15.jpg");
		Matrix test = getMatrixFromFile("C://Users//STEVEN-PC//Desktop//centerimgs//subject05.normal");

		
		Matrix meanSub = test.minus(meanVector);
		Matrix testomega = covEigVec.transpose().times(meanSub);
		
		Matrix mat = reconstructFace(testomega, covEigVec);
		saveEigenfaces(mat, "C://Users//STEVEN-PC//Desktop//save");

		double a1 = getMinError(test, meanVector, covEigVec, omega); 
			
		System.out.println(a1);
		
//		System.out.println(a1);
		
		saveEigenfaces(covEigVec, "C://Users//STEVEN-PC//Desktop//eigenfaces//pic");
	}
	
	//reconstruct face given weights
	public static Matrix reconstructFace(Matrix weights, Matrix eigenVec){
		Matrix sum = new Matrix(eigenVec.getRowDimension(), 1);
		for (int i = 0; i < eigenVec.getColumnDimension(); i++){
			Matrix curr = getColumn(eigenVec, i);
			double weight = weights.get(i, 0);
			sum = sum.plus(curr.times(weight));
		}
		
		return sum;
	}
	
	// get minimum error given Matrix m of test image
	public static double getMinError(Matrix m, Matrix mean, Matrix eigVec, Matrix omega){
		Matrix meanSub = m.minus(mean);
		Matrix testomega = normalizeColumns(eigVec.transpose().times(meanSub));
		//compare with columns of omega with euclidean distance to find minimum
		int numCols = omega.getColumnDimension();
		double minimum = Double.MAX_VALUE;
		for (int i = 0 ; i < numCols; i++){
			Matrix col = normalizeColumns(getColumn(omega, i));
			Matrix diff = testomega.minus(col);
			double norm = diff.normF();
			if (norm < minimum)
				minimum = norm;
		}
		
		return minimum;
	}
	
	// get column from matrix given index
	public static Matrix getColumn(Matrix m, int i){
		int numRows = m.getRowDimension();
		double[] data = new double[numRows];
		for (int a = 0 ; a < numRows; a++)
			data[a] = m.get(a, i);
		Matrix mat = new Matrix(data, numRows);
		return mat;
	}
	
	// omega = U' * (L - mean)
	public static Matrix getWeights(Matrix eigVec, Matrix diffMat){
		return eigVec.transpose().times(diffMat);
	}
	
	// retrieve M * N vector from image file
	public static Matrix getMatrixFromFile(String dir){
		File f = new File(dir);
		ImagePlus img = new ImagePlus(f.toString());
		BufferedImage bImg = img.getBufferedImage();
		double[] db = new double[DIMENSION];
		int index = 0;
		for (int i = 0; i < bImg.getHeight(); i++){
			for (int j = 0; j < bImg.getWidth(); j++){
				db[index] = bImg.getRGB(j, i)&0xFF;
				index++;
			}
		}
		Matrix m = new Matrix(db, DIMENSION);
		return m;
	}
	
	//AAT = covariance matrix; ATA has same eigenvalues and Avi are eigenvectors of AAT; U is normalized
	public static Matrix getEigenface(Matrix differenceMatrix){
		//compute ATA
		Matrix med = differenceMatrix.transpose().times(differenceMatrix);
		EigenvalueDecomposition decomp = med.eig();
		Matrix aEigVal = decomp.getD();
		Matrix aEigVec = decomp.getV();
		Matrix covEigVec = differenceMatrix.times(aEigVec); //u = Av
		Matrix normCovEigVec = normalizeColumns(covEigVec);
		normCovEigVec = flipMatrix(normCovEigVec);
		return normCovEigVec;
	}
	
	//normalizes columns
	public static Matrix normalizeColumns(Matrix m){
		int numRows = m.getRowDimension();
		int numCols = m.getColumnDimension();
		
		Matrix normalizedM = new Matrix(numRows, numCols);
		for (int i = 0; i < numCols; i++){
			double sum = 0;
			for (int j = 0; j < numRows; j++)
				sum += Math.pow(m.get(j, i),2);
			sum = Math.sqrt(sum);
			for (int j = 0; j < numRows; j++)
				normalizedM.set(j, i, m.get(j, i)/sum);
		}
		
		return normalizedM;
	}
	
	// flip columns to opposite side 
	public static Matrix flipMatrix(Matrix m){
		int cols = m.getColumnDimension();
		int rows = m.getRowDimension();
		
		Matrix a = new Matrix(rows, cols);
		for (int i = 0; i < cols; i++){
			for (int j = 0; j < rows; j++){
				a.set(j, cols - (i+1), m.get(j, i));
			}
		}
		
		return a;
	}
	
	//save eigenfaces into image of same dimensions
	public static void saveEigenfaces(Matrix covEigVec, String dir) throws Exception{
		double[][] data = covEigVec.transpose().getArray();
		
		// matrix is reverse sorted
		for (int i = covEigVec.getColumnDimension()-1; i >= 0; i--){
			
			double[] row = data[i];
			List list = Arrays.asList(ArrayUtils.toObject(row));
			double min = Collections.min(list);
			double max = Collections.max(list);
			
			for (int j = 0; j < row.length; j++)
				row[j] = getScaledValue(row[j], min, max);

			BufferedImage newimg = new BufferedImage(320, 243, BufferedImage.TYPE_INT_RGB);
			populateImage(newimg, row);
			
			File outputFile = new File(dir + i + ".jpg");
			ImageIO.write(newimg, "jpg", outputFile);
		}
	}

	//populate bufferedimage using double[] with pixel values
	public static void populateImage(BufferedImage img, double[] data){
		int a = 0;
		for (int i = 0; i < img.getHeight(); i++){
			for (int j = 0 ; j < img.getWidth(); j++){
				int currRGB = (int) data[a];
				Color col = new Color(currRGB, currRGB, currRGB);
				int RGB = col.getRGB();
				img.setRGB(j, i, RGB);
				a++;
			}
		}
	}
	
	//get difference matrix (mean - original)
	public static Matrix getDifferenceMatrix(Matrix meanVector, Matrix dataMatrix){
		double[][] meanArrayVector2D = meanVector.getArrayCopy();
		double[] meanArrayVector = new double[meanVector.getRowDimension()];
		double[][] meanMatrixExpanded = new double[dataMatrix.getColumnDimension()][dataMatrix.getRowDimension()];
		
		for (int i = 0; i < meanVector.getRowDimension(); i++)
			meanArrayVector[i] = meanArrayVector2D[i][0];
		
		for (int i = 0; i < dataMatrix.getColumnDimension(); i++)
			meanMatrixExpanded[i] = meanArrayVector;
		
		Matrix totalMean = new Matrix(meanMatrixExpanded);
		totalMean = totalMean.transpose();
				
		Matrix differenceMatrix = dataMatrix.minus(totalMean);
		return differenceMatrix;
	}

	//convert pixel to 0-255 range based on range in eigenvector
	public static double getScaledValue(double pixelValue, double min, double max){
		return 255 * (pixelValue - min)/(max-min);
	}
	
	//get Matrix with data of pixel values
	public static Matrix getDataMatrix(File[] files){
		int len = files.length;
//		int len = 154;
		double[][] data = new double[DIMENSION][len];
		for (int i = 0; i < len; i++){
			File file = files[i];
			ImagePlus img = new ImagePlus(file.toString());
			BufferedImage bImg = img.getBufferedImage();
			int index = 0;
			for (int j = 0; j < bImg.getHeight(); j++){
				for (int k = 0; k < bImg.getWidth(); k++){
					data[index][i] = bImg.getRGB(k, j)&0xFF;
					index++;
				}
			}
		}
		
		Matrix matr = new Matrix(data);
		return matr;
	}
	
	//find vector with mean pixel values
	public static Matrix getMeanVector(Matrix total){
		double[] meanVector = new double[total.getRowDimension()];
		
		for (int i = 0; i < total.getColumnDimension(); i++)
			for (int j = 0; j < total.getRowDimension(); j++)
				meanVector[j] += total.get(j, i);
		
		for (int i = 0; i < total.getRowDimension(); i++)
			meanVector[i] /= total.getColumnDimension();	
		
		Matrix mat = new Matrix(meanVector, total.getRowDimension());  
		return mat;
	}
	
	//get files from directory
	public static File[] getFiles(String baseDir){
		File[] files = new File(baseDir).listFiles();
		return files;
	}
}
