package evaluation;

import java.io.File;

import survivalFactorizationEM.SurvivalFactorizationEM_Model;
import utils.MatrixUtilities;

public class OriginalMatrixPrinter {
	public static void main(String[] args) throws Exception {
		final String modelFolder = args[0];
		final String outputFolder = args[1];

		final File[] files = new File(modelFolder).listFiles();

		for (final File f : files)
			if (f.getName().endsWith(".model")) {
				System.out.print("Loading model" + f.getName() + "...");
				final SurvivalFactorizationEM_Model model = SurvivalFactorizationEM_Model
						.readFromFile(f.getAbsolutePath());
				System.out.println("... done!");

				System.out.println("Matrix A");
				MatrixUtilities.print(model.A, outputFolder + File.separator
						+ f.getName() + ".A");

				System.out.println("Matrix S");
				MatrixUtilities.print(model.S, outputFolder + File.separator
						+ f.getName() + ".S");

				if (model.Phi != null) {
					System.out.println("Matrix Phi");
					MatrixUtilities.print(model.Phi, outputFolder
							+ File.separator + f.getName() + ".phi");
				}

				System.out.println("\r\n... done!");
				System.out.println("\r\n");
			}
	}
}