import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Date;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.*;
import java.util.ArrayList;

public class predictBatchPEPSTATS
{
  private static final String EMBOSS_VERSION = "EMBOSS-6.5.7";
  private static final String LIBSVM_VERSION = "libsvm-3.16";
  
  public static void main(String[] args, Vector<String> ecnums, long time, String ROOTPATH, String fastaFile, String tempDir)
    throws IOException, InterruptedException
  {
    String method = "pepstats";
    
    // Parallelize EC processing
    ParallelExecutor executor = ParallelExecutor.getInstance();
    List<Callable<Void>> ecTasks = new ArrayList<>();
    
    for (int i = 0; i < ecnums.size(); i++)
    {
      final int ecIndex = i;
      final String ecnum = ecnums.get(ecIndex);
      
      ecTasks.add(() -> {
        try {
          Date d1 = new Date();
          
          String path = ROOTPATH + File.separator + ecnum + File.separator + method;
          String modelfile = path + File.separator + "model.svm";
          String rangefile = path + File.separator + "rangefile";
          String testpath = tempDir + File.separator + "testResult" + File.separator + time + File.separator + ecnum + File.separator + method;
          
          String batchSVM = testpath + "/test.svm";
          String batchVect = testpath + "/test.vec";
          File workdir = new File(testpath);
          workdir.mkdirs();
          
          fasta2Pepstats_noscale fas = new fasta2Pepstats_noscale();
          
          String predFile = testpath + File.separator + ecnum + ".preds";
          String confFile = testpath + File.separator + ecnum + ".confs";
          
          String posPredFile = path + File.separator + "ppreds.txt";
          String negPredFile = path + File.separator + "npreds.txt";
          String[] cmdArray = new String[7];
          
          cmdArray[0] = (ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/" + EMBOSS_VERSION + "/emboss/pepstats");
          
          cmdArray[1] = "-sequence";
          
          cmdArray[2] = fastaFile;
          
          cmdArray[3] = "-outfile";
          
          cmdArray[4] = (workdir + "/temp.out");
          
          cmdArray[5] = "-warning";
          
          cmdArray[6] = "FALSE";
          
          String cmd = ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/" + EMBOSS_VERSION + "/emboss/pepstats -sequence " + fastaFile + " -outfile " + workdir + "/temp.out -warning FALSE";
          ProcessBuilder pb = new ProcessBuilder(cmd.split(" "));
          
          Process process = pb.start();
          try
          {
            process.waitFor();
            int exitVal = process.exitValue();
          }
          catch (InterruptedException e)
          {
            System.out.print("pepstats is not working!");
            e.printStackTrace();
          }
          fasta2Pepstats_noscale.parse_pepstats(cmdArray[4], ecnum, "1", time, ROOTPATH, tempDir);
          File deletefile = new File(cmdArray[4]);
          deletefile.delete();
          
          // Run svm-scale directly without shell script to avoid race conditions
          String scaleCmd = ROOTPATH.substring(0, ROOTPATH.length() - 3) + "/" + LIBSVM_VERSION + "/svm-scale";
          String rangeFileArg = ROOTPATH + "/" + ecnum + "/pepstats/rangefile";
          
          pb = new ProcessBuilder(scaleCmd, "-r", rangeFileArg, batchSVM);
          pb.redirectOutput(new File(batchVect));
          pb.redirectError(ProcessBuilder.Redirect.DISCARD);
          
          process = pb.start();
          try
          {
            process.waitFor();
            int k = process.exitValue();
          }
          catch (InterruptedException e)
          {
            System.out.print("svm-scale is not working!");
            e.printStackTrace();
          }
          
          // Use in-JVM SVM classifier instead of external process
          SVMLightClassifier.classify(batchVect, modelfile, predFile);
          
          utils u = new utils();
          utils.calculateConfidence(posPredFile, negPredFile, predFile, confFile);
        } catch (Exception e) {
          System.err.println("Error processing EC " + ecnum + " in PEPSTATS: " + e.getMessage());
          e.printStackTrace();
        }
        return null;
      });
    }
    
    // Execute all EC tasks in parallel
    try {
      executor.executeECClassLevel(ecTasks);
    } catch (ExecutionException e) {
      System.err.println("Error executing parallel PEPSTATS EC processing: " + e.getMessage());
      throw new IOException("Parallel PEPSTATS execution failed", e);
    }
  }
}