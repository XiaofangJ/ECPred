import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

public class runEC {

	private static final String[] METHODS = {"blast", "spmap", "pepstats"};
	private static final double NON_ENZYME_THRESHOLD = 0.4;

	public HashMap<String, Vector<Vector<String>>> predictions(String[] args, String ROOTPATH, Vector<String> ecnums, long time, HashMap<String, Vector<Vector<String>>> predictions, List<String> idlist, String fastaFile, String tempDir, String method) throws IOException, InterruptedException{	
		
		// Run predictions based on method
		runPredictionMethods(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir, method);
		
		// Load thresholds
		HashMap<String, Double> thresholds = loadThresholds(ROOTPATH);
		
		// Process results for each EC number
		for (int i = 0; i < ecnums.size(); i++) {
			String testDir = tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + "preds";
			createTestDirectory(testDir);
			
			Vector<String> combined = loadAndCombinePredictions(method, tempDir, time, ecnums.get(i), ROOTPATH, idlist);
			writePredictionFile(combined, tempDir, time, ecnums.get(i), idlist);
		}
		
		// Process main class or subclass predictions
		if (ecnums.contains("1.-.-.-")) {
			processMainClassPredictions(ecnums, tempDir, time, ROOTPATH, idlist, thresholds, predictions);
		} else {
			processSubclassPredictions(ecnums, tempDir, time, thresholds, predictions, idlist);
		}
		
		return predictions;
	}

	private void runPredictionMethods(String[] args, Vector<String> ecnums, long time, String ROOTPATH, List<String> idlist, String fastaFile, String tempDir, String method) throws IOException, InterruptedException {
		if (method.equals("spmap")) {
			predictBatchSPMAP.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
		} else if (method.equals("blast")) {
			predictBatchBLAST.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
		} else if (method.equals("pepstats")) {
			predictBatchPEPSTATS.main(args, ecnums, time, ROOTPATH, fastaFile, tempDir);
		} else if (method.equals("weighted")) {
			// Run all three methods for weighted approach
			predictBatchSPMAP.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
			predictBatchBLAST.main(args, ecnums, time, ROOTPATH, idlist, fastaFile, tempDir);
			predictBatchPEPSTATS.main(args, ecnums, time, ROOTPATH, fastaFile, tempDir);
		}
	}

	private HashMap<String, Double> loadThresholds(String ROOTPATH) throws IOException {
		HashMap<String, Double> thresholds = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/subclasses/thresholds.txt"));
		String line;
		while ((line = br.readLine()) != null) {
			StringTokenizer st1 = new StringTokenizer(line, "\t");
			String ECclass = st1.nextToken();
			String thres = st1.nextToken();
			Double threshold = Double.valueOf(Double.parseDouble(thres));
			thresholds.put(ECclass, threshold);
		}
		br.close();
		return thresholds;
	}

	private void createTestDirectory(String testDir) {
		File folder1 = new File(testDir);
		folder1.mkdirs();
		if (folder1.listFiles() != null) {
			for (final File fileEntry : folder1.listFiles()) {
				fileEntry.delete();
			}
		}
	}

	private Vector<String> loadAndCombinePredictions(String method, String tempDir, long time, String ecnum, String ROOTPATH, List<String> idlist) throws IOException {
		Vector<String> combined = new Vector<>();
		
		if (method.equals("weighted")) {
			return loadWeightedPredictions(tempDir, time, ecnum, ROOTPATH, idlist);
		} else {
			return loadSingleMethodPredictions(method, tempDir, time, ecnum);
		}
	}

	private Vector<String> loadSingleMethodPredictions(String method, String tempDir, long time, String ecnum) throws IOException {
		Vector<String> combined = new Vector<>();
		String predFile = tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + method + File.separator + ecnum + ".confs";
		List<String> preds = Files.readAllLines(Paths.get(predFile));
		
		for (String pred : preds) {
			Double comb = Double.valueOf(Double.parseDouble(pred));
			combined.add(String.valueOf(comb));
		}
		return combined;
	}

	private Vector<String> loadWeightedPredictions(String tempDir, long time, String ecnum, String ROOTPATH, List<String> idlist) throws IOException {
		Vector<String> combined = new Vector<>();
		
		List<String> spreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + "spmap" + File.separator + ecnum + ".confs"));
		List<String> bpreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + "blast" + File.separator + ecnum + ".confs"));
		List<String> ppreds = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + "pepstats" + File.separator + ecnum + ".confs"));
		List<String> weights = Files.readAllLines(Paths.get(ROOTPATH + File.separator + ecnum + File.separator + "weights.txt"));

		for (int j = 0; j < spreds.size(); j++) {
			Double comb = Double.valueOf(Double.parseDouble(spreds.get(j)) * Double.parseDouble(weights.get(0)) + 
										  Double.parseDouble(bpreds.get(j)) * Double.parseDouble(weights.get(1)) + 
										  Double.parseDouble(ppreds.get(j)) * Double.parseDouble(weights.get(2)));
			combined.add(String.valueOf(comb));
		}

		// Write detailed predictions for weighted method
		String testDir = tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + "preds";
		for (int j = 0; j < idlist.size(); j++) {
			BufferedWriter final_file1 = new BufferedWriter(new FileWriter(testDir + File.separator + idlist.get(j) + ".preds", true));
			final_file1.write(ecnum + "\t" + spreds.get(j) + "\t" + bpreds.get(j) + "\t" + ppreds.get(j) + "\t" + combined.get(j) + "\n");
			final_file1.close();
		}
		
		return combined;
	}

	private void writePredictionFile(Vector<String> combined, String tempDir, long time, String ecnum, List<String> idlist) throws IOException {
		DecimalFormat df = new DecimalFormat();
		df.setMaximumFractionDigits(2);
		BufferedWriter final_file = new BufferedWriter(new FileWriter(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + ecnum + "_preds.txt", false));
		for (int j = 0; j < idlist.size(); j++) {
			final_file.write(df.format(Double.parseDouble(combined.get(j))) + "\n");
		}
		final_file.close();
	}

	private void processMainClassPredictions(Vector<String> ecnums, String tempDir, long time, String ROOTPATH, List<String> idlist, HashMap<String, Double> thresholds, HashMap<String, Vector<Vector<String>>> predictions) throws IOException {
		Vector<Vector<String>> allPreds = new Vector<>();
		
		// Load all predictions
		for (int i = 0; i < ecnums.size(); i++) {
			Vector<String> pred = new Vector<>();
			BufferedReader br = new BufferedReader(new FileReader(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + ecnums.get(i) + "_preds.txt"));
			String line;
			while ((line = br.readLine()) != null) {
				pred.add(line);
			}
			br.close();
			allPreds.add(pred);
		}

		// Process each protein
		for (int i = 0; i < idlist.size(); i++) {
			Vector<Vector<String>> predswithScore = new Vector<>();
			Vector<String> preds = new Vector<>();
			double maxPred = 0.0;
			String mainClass = null;
			
			// Find the class with highest prediction
			for (int j = 0; j < allPreds.size(); j++) {
				double currentPred = Double.parseDouble(allPreds.get(j).get(i));
				if (currentPred > maxPred) {
					maxPred = currentPred;
					mainClass = String.valueOf((j + 1)) + ".-.-.-";
				}
			}

			// Classify based on thresholds
			if (maxPred < NON_ENZYME_THRESHOLD) {
				preds.add("non");
				preds.add(String.valueOf(1.0 - maxPred));
			} else if (maxPred >= thresholds.get(mainClass)) {
				preds.add(mainClass);
				preds.add(String.valueOf(maxPred));
			} else {
				preds.add("nop");
				preds.add("0");
			}
			
			predswithScore.add(preds);
			predictions.put(idlist.get(i), predswithScore);
		}
	}

	private void processSubclassPredictions(Vector<String> ecnums, String tempDir, long time, HashMap<String, Double> thresholds, HashMap<String, Vector<Vector<String>>> predictions, List<String> idlist) throws IOException {
		double maxPred = 0.0;
		Vector<String> preds = new Vector<>();
		
		// Find best subclass prediction
		for (int i = 0; i < ecnums.size(); i++) {
			List<String> pred = Files.readAllLines(Paths.get(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnums.get(i) + File.separator + ecnums.get(i) + "_preds.txt"));
			String predClass = ecnums.get(i);
			double currentPred = Double.parseDouble(pred.get(0));
			
			if (currentPred >= thresholds.get(predClass) && currentPred > maxPred) {
				maxPred = currentPred;
				preds = new Vector<>();
				preds.add(predClass);
				preds.add(pred.get(0));
			}
		}
		
		if (preds.size() == 0) {
			preds.add("nop");
			preds.add("0");
		}

		// Add prediction to existing protein's prediction vector
		if (predictions.containsKey(idlist.get(0))) {
			predictions.get(idlist.get(0)).add(preds);
		} else {
			Vector<Vector<String>> predswithScore = new Vector<>();
			predswithScore.add(preds);
			predictions.put(idlist.get(0), predswithScore);
		}
	}

}