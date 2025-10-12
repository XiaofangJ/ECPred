import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;

public class predictBatchSPMAP {

	public static void main(String[] args, Vector<String> ecnums, long time, String ROOTPATH, List<String> test_ids, String fastaFile, String tempDir)
		    throws IOException, InterruptedException
		  {
		    int sigTh = -15;
		    int subseqlen = 5;
		    String method = "spmap";
		    
		    // Parallelize EC processing
		    ParallelExecutor executor = ParallelExecutor.getInstance();
		    List<Callable<Void>> ecTasks = new ArrayList<>();
		    
		    for (int i = 0; i < ecnums.size(); i++)
		    {
		      final int ecIndex = i;
		      final String ecnum = ecnums.get(ecIndex);
		      
		      ecTasks.add(() -> {
		        try {
		          File workdir = new File(tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum);
		          workdir.mkdirs();
		          String path = ROOTPATH + File.separator + ecnum + File.separator + method;
		          String testpath = workdir + File.separator + method;
		          String modelfile = path + File.separator + "model.svm";
		          
		          workdir = new File(testpath);
		          workdir.mkdirs();
		          
		          String predFile = testpath + File.separator + ecnum + ".preds";
		          String confFile = testpath + File.separator + ecnum + ".confs";
		          String batchVect = testpath + File.separator + "test.vec";
		          
		          String posPredFile = path + File.separator + "ppreds.txt";
		          String negPredFile = path + File.separator + "npreds.txt";
		          
		          seq2vectPSSMtest.calculateVectors(sigTh, subseqlen, ecnum, test_ids, fastaFile, time, ROOTPATH, tempDir);
		          
		          Path batchVectPath = Paths.get(batchVect);
		          if (Files.notExists(batchVectPath) || Files.size(batchVectPath) == 0) {
		            System.err.println("Warning: Skipping SPMAP classification for EC " + ecnum + " due to empty feature vector file.");
		            return null;
		          }
		          
		          // Use in-JVM SVM classifier instead of external process
		          SVMLightClassifier.classify(batchVect, modelfile, predFile);
		          
		          Path predictionPath = Paths.get(predFile);
		          if (Files.exists(predictionPath) && Files.size(predictionPath) > 0) {
		            utils.calculateConfidence(posPredFile, negPredFile, predFile, confFile);
		          } else {
		            System.err.println("Warning: Missing prediction output for EC " + ecnum + ", skipping confidence calculation.");
		          }
		        } catch (Exception e) {
		          System.err.println("Error processing EC " + ecnum + " in SPMAP: " + e.getMessage());
		          e.printStackTrace();
		        }
		        return null;
		      });
		    }
		    
		    // Execute all EC tasks in parallel
		    try {
		      executor.executeECClassLevel(ecTasks);
		    } catch (ExecutionException e) {
		      System.err.println("Error executing parallel SPMAP EC processing: " + e.getMessage());
		      throw new IOException("Parallel SPMAP execution failed", e);
		    }
	}
}