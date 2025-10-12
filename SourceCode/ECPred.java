/*ECPred: a tool for the prediction of enzymatic properties of protein sequences based on the EC Nomenclature
    Copyright (C) 2018 CanSyL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.concurrent.*;
import java.util.concurrent.ConcurrentHashMap;

public class ECPred {

	
		
	public static void main(String[] args) throws IOException, InterruptedException {
		if (args.length < 4)
	    {
	      System.out.println("Missing Argument(s)!");
	      System.out.println("Sample run: java -jar ECPred.jar method inputFile libraryDir tempDir outputFile");
	      System.out.println("method argument can be one of the followings: blast, spmap, pepstats, weighted");
	      System.out.println("outputFile is optinal. If you don't specify the output file name the results will be printed to standard output.");
	      System.exit(0);
	    }
		Date d1 = new Date();
		long time = d1.getTime();
		long t0 = System.currentTimeMillis();
	    String output = "";String tempDir = "";
	    Vector<String> ecnums = new Vector<>();
	    ecnums.add("1.-.-.-");ecnums.add("2.-.-.-");ecnums.add("3.-.-.-");
	    ecnums.add("4.-.-.-");ecnums.add("5.-.-.-");ecnums.add("6.-.-.-");
	    runEC runECobj = new runEC();
	    
	    String method = args[0];
	    if ((!method.equals("blast")) && (!method.equals("spmap")) && (!method.equals("pepstats")) && (!method.equals("weighted")))
	    {
	      System.out.println("method argument must be one of the followings: blast, spmap, pepstats, weighted");
	      System.exit(0);
	    }
	    String fastaFile = args[1];
	    String ROOTPATH = args[2] + "lib/EC";
	    String file_namepart1 = "";
	    String file_name;
	    if (fastaFile.contains("/"))
	    {
	      StringTokenizer st1 = new StringTokenizer(fastaFile, "/");
	      while (st1.hasMoreTokens()) {
	        file_namepart1 = st1.nextToken();
	      }
	      st1 = new StringTokenizer(file_namepart1, ".");
	      file_name = st1.nextToken();
	    }
	    else
	    {
	      StringTokenizer st1 = new StringTokenizer(fastaFile, ".");
	      file_name = st1.nextToken();
	    }
	    
		// Robust argument parsing for optional CPU argument
		int userCpuCount = Runtime.getRuntime().availableProcessors();
		System.out.println("[TIMER] Program started at: " + new Date());
	    if (args.length == 4) {
	        output = "stdout";
	        tempDir = args[3];
	    } else if (args.length == 5) {
	        output = args[4];
	        tempDir = args[3];
	    } else if (args.length == 6) {
	        output = args[4];
	        tempDir = args[3];
	        try {
	            userCpuCount = Integer.parseInt(args[5]);
	            if (userCpuCount < 1) userCpuCount = Runtime.getRuntime().availableProcessors();
	        } catch (Exception e) {
	            userCpuCount = Runtime.getRuntime().availableProcessors();
	        }
	    }
	    ParallelExecutor.setUserCpuCount(userCpuCount);
		
		String dateandtime = new SimpleDateFormat("yyyyMMdd-HHmmss").format(new Date());

		long t2 = System.currentTimeMillis();
		HashMap<String, String> protID = checkFasta(fastaFile);
		List<String> idlist = createFasta(fastaFile);
		System.out.println("[TIMER] FASTA parsing: " + (System.currentTimeMillis() - t2) + " ms");


		ConcurrentHashMap<String, Vector<Vector<String>>> predictions = new ConcurrentHashMap<>();
		System.out.println("Main classes of input proteins are being predicted ...");
		
		final String finalTempDir = tempDir;
		final long finalTime = time;
		final String finalFastaFile = fastaFile;
		final String finalROOTPATH = ROOTPATH;
		final String finalMethod = method;
		final String[] finalArgs = args;
		
		createFasta(idlist, fastaFile, "test.fasta", tempDir+File.separator + "testResult" +File.separator + time);
		String newfasta =  tempDir + File.separator + "testResult" + File.separator + time + File.separator + "test.fasta"; 
		
	// First pass: predict main classes for all proteins
	int blastThreads = Math.max(1, userCpuCount / 6);
	
	long t4 = System.currentTimeMillis();
	HashMap<String, Vector<Vector<String>>> mainPredictions = runECobj.predictions(args, ROOTPATH, ecnums, time, new HashMap<>(), idlist, newfasta, tempDir, method, blastThreads);
	predictions.putAll(mainPredictions);
		
	// Second pass: predict subclasses for each protein in parallel
	// Use ParallelExecutor for consistent thread pool management
	ParallelExecutor executor = ParallelExecutor.getInstance();
	List<CompletableFuture<Void>> proteinFutures = new ArrayList<>();

	// Submit protein tasks using ParallelExecutor
	long t5 = System.currentTimeMillis();
	for (Map.Entry<String, Vector<Vector<String>>> entry : predictions.entrySet()) {
		if (entry.getValue().get(0).get(0).equals("non") || entry.getValue().get(0).get(0).equals("nop"))
			continue;

		final String proteinKey = entry.getKey();
		final String proteinTempDir = tempDir + File.separator + "testResult" + File.separator + finalTime + "_" + proteinKey.hashCode();
		
		CompletableFuture<Void> proteinFuture = executor.submitProteinTask(() -> {
			try {
				synchronized(System.out) {
					if(protID.get(proteinKey).length()>81)
						System.out.println("Subclasses of "+protID.get(proteinKey).substring(1,81) + " are being predicted ...");
					else
						System.out.println("Subclasses of "+protID.get(proteinKey).substring(1,protID.get(proteinKey).length()) + " are being predicted ...");
				}
				
				for(int i = 1 ; i<4; i++){
					List<String> ecList = Files.readAllLines(Paths.get(finalROOTPATH.substring(0, finalROOTPATH.length()-3)+"/subclasses/"+ predictions.get(proteinKey).get(i-1).get(0) + ".txt"));
					if(ecList.size()==0){
						Vector<String> preds = new Vector<>();
						preds.add("nop");
						preds.add("");
						synchronized(predictions.get(proteinKey)) {
							predictions.get(proteinKey).add(preds);
						}
						break;
					}

					Vector<String> ecnumsLocal = new Vector<>() ;
					ecnumsLocal.addAll(ecList);
					List<String> idlistLocal = new Vector<String>();
					idlistLocal.add(proteinKey);
					
					synchronized(ECPred.class) {
						createFasta(idlistLocal, finalFastaFile, "test.fasta", proteinTempDir);
					}
					final String newfastaLocal = proteinTempDir + File.separator + "test.fasta";
					
					HashMap<String, Vector<Vector<String>>> subPredictions = runECobj.predictions(finalArgs, finalROOTPATH, ecnumsLocal, finalTime, new HashMap<>(), idlistLocal, newfastaLocal, proteinTempDir, finalMethod, 1);
					
					if(subPredictions.containsKey(proteinKey)) {
						Vector<Vector<String>> subPreds = subPredictions.get(proteinKey);
						if(subPreds.size() > 0) {
							synchronized(predictions.get(proteinKey)) {
								predictions.get(proteinKey).addAll(subPreds);
							}
						}
					}
					
					if(predictions.get(proteinKey).size() > i && predictions.get(proteinKey).get(i).get(0).equals("nop"))
						break;
				}
			} catch (Exception e) {
				System.err.println("Error predicting subclasses for protein " + proteinKey + ": " + e.getMessage());
				e.printStackTrace();
			}
			return null;
		});

		proteinFutures.add(proteinFuture);
	}
	
	// Wait for all protein predictions to complete
	try {
		CompletableFuture.allOf(proteinFutures.toArray(new CompletableFuture[0])).join();
	} catch (Exception e) {
		System.err.println("Error in parallel protein processing: " + e.getMessage());
		e.printStackTrace();
	}
	System.out.println("[TIMER] Subclass prediction: " + (System.currentTimeMillis() - t5) + " ms");
	// Now collect results for output, ensuring all subclass predictions are present
Map<String, List<String>> finalResults = new ConcurrentHashMap<>();
for (String id : idlist) {
    if (!predictions.containsKey(id)) continue;
    Vector<Vector<String>> predVec = predictions.get(id);
    StringBuilder sb = new StringBuilder();
    if (protID.get(id).length() > 81)
        sb.append(protID.get(id).substring(1, 81));
    else
        sb.append(protID.get(id).substring(1, protID.get(id).length()));
	if (predVec.get(0).get(0).equals("non")) {
		double conf = 1.0 - Double.parseDouble(predVec.get(0).get(1));
		sb.append("\tnon Enzyme\t").append(String.format("%.2f", conf));
    } else if (predVec.get(0).get(0).equals("nop")) {
        sb.append("\tno Prediction");
    } else {
        // Find the most specific (deepest) prediction that is not "nop"
        int bestIdx = 0;
        for (int i = 1; i < predVec.size(); i++) {
            if (!predVec.get(i).get(0).equals("nop") && predVec.get(i).get(0).contains(".")) {
                bestIdx = i;
            } else {
                break;
            }
        }
        sb.append("\t").append(predVec.get(bestIdx).get(0)).append("\t").append(predVec.get(bestIdx).get(1));
    }
    finalResults.put(id, List.of(sb.toString()));
}

	long t6 = System.currentTimeMillis();
	if (!output.equals("stdout")) {
		PrintWriter predFile = new PrintWriter(output, "UTF-8");
		predFile.println("Protein ID\tEC Number\tConfidence Score(max 1.0)");
		for (String id : idlist) {
			if (finalResults.containsKey(id)) {
				predFile.println(finalResults.get(id).get(0));
			}
		}
		predFile.close();
	} else {
		System.out.println("Protein ID\tEC Number\tConfidence Score(max 1.0)");
		for (String id : idlist) {
			if (finalResults.containsKey(id)) {
				System.out.println(finalResults.get(id).get(0));
			}
		}
	}
	System.out.println("[TIMER] Output writing: " + (System.currentTimeMillis() - t6) + " ms");
			
	    		 
			Date d2 = new Date();
			long diff = d2.getTime() - d1.getTime();
			long diffSeconds = diff / 1000 % 60;
			long diffMinutes = diff / (60 * 1000) % 60;
			long diffHours = diff / (60 * 60 * 1000) % 24;
			long diffDays = diff / (24 * 60 * 60 * 1000);

			System.out.println("[TIMER] Total runtime: " + (System.currentTimeMillis() - t0) + " ms");
			if(diffMinutes==0 && diffHours ==0 && diffDays==0)
				System.out.println("--- Proteins are predicted in "+diffSeconds + " seconds ---");
			else if(diffHours ==0 && diffDays==0)
				System.out.println("--- Proteins are predicted in "+diffMinutes + " minutes "+diffSeconds + " seconds ---");
			else if(diffDays==0)
				System.out.println("--- Proteins are predicted in "+diffHours + " hours "+diffMinutes + " minutes "+diffSeconds + " seconds ---");
				else
					System.out.print("--- Proteins are predicted in "+diffDays + " days "+diffHours + " hours "+diffMinutes + " minutes "+diffSeconds + " seconds ---");
				
				// Shutdown parallel executor
				ParallelExecutor.getInstance().shutdown();	}
	
	public static boolean deleteDirectory(File directory) {
	    if(directory.exists()){
	        File[] files = directory.listFiles();
	        if(null!=files){
	            for(int i=0; i<files.length; i++) {
	                if(files[i].isDirectory()) {
	                    deleteDirectory(files[i]);
	                }
	                else {
	                    files[i].delete();
	                }
	            }
	        }
	    }
	    return(directory.delete());
	}
	
	public static HashMap<String, String> checkFasta( String fasta) throws IOException{
		HashMap<String, String> protID_Full = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(fasta));	
		List<String> idlist = new ArrayList<String>();
		String orjinal_line,line = br.readLine(), fastaString = "";
		String prot_id;
		 while(line != null ){
			 fastaString = "";
			 if(line.startsWith(">")==false){
					System.out.println("Wrong input! Sequences should start with \">\" character.");
					System.exit(0);
				}
			 orjinal_line = line;
			 line = line.replaceAll("[^\\dA-Za-z ]", "").replaceAll("\\s+", "+");
			   line = line.trim().replaceAll(" +", "");
			   if(line.length()>80)
					prot_id = line.substring(1, 81);
			   else
				   prot_id  = line.substring(1, line.length());
			 protID_Full.put(prot_id, orjinal_line);
			line = br.readLine();
		    while(line.startsWith(">")==false){ // Read sequence lines
		    	 fastaString+=line;
				line = br.readLine();
				if(line==null)
					 break;
			   }
		    fastaString = fastaString.replaceAll("\\s+","");
		    Pattern p = Pattern.compile("[^A-Z ]", Pattern.CASE_INSENSITIVE);
		      Matcher m = p.matcher(fastaString);
		      boolean flag =false;
		      while (m.find()) {
		         System.out.println("Fasta sequence contains special character at position "  + m.start() + ": " + fastaString.charAt(m.start()) + " Your fasta: "+fastaString);
		         flag = true;
		      }
		      if(flag== true){
			      System.exit(0);
		      }
		}	
		 br.close();
		return protID_Full;
		
	}

	
	public static List<String> createFasta( String fasta) throws IOException{
		String prot_id = "";
		BufferedReader br = new BufferedReader(new FileReader(fasta));	
		List<String> idlist = new ArrayList<String>();
		String line = br.readLine(), fastaString = "";
		
		 while(line != null ){
				 fastaString = "";
			
				   line = line.replaceAll("[^\\dA-Za-z ]", "").replaceAll("\\s+", "+");
				   line = line.trim().replaceAll(" +", "");
				  if(line.length()>80)
						prot_id = line.substring(1, 81);
				   else
					   prot_id  = line.substring(1, line.length());
				line = br.readLine();
			    while(line.startsWith(">")==false){//until > add sequence 
			    	 fastaString+=line;
					line = br.readLine();
					if(line==null)
						 break;
				   }
			    if(fastaString.length()>40)
					idlist.add(prot_id);
			 	
			}	
			 br.close();	
		 return idlist;
	}

	public static void createFasta( List<String> lst_training_ids,String fasta, String outf, String file) throws IOException{
		File workdir = new File(file);
		workdir.mkdirs();
		PrintWriter final_file = new PrintWriter(file+"/"+outf, "UTF-8");
		String prot_id = "";
		BufferedReader br = new BufferedReader(new FileReader(fasta));	
		HashMap<String, String> fastaArray = new HashMap<>();
		Vector<String> fastaWrite = new Vector<>();
		StringTokenizer st1;
		String line = br.readLine(), fastaString = "";
		int i=0;
		 while(line != null ){
			 fastaString = "";
			   line = line.replaceAll("[^\\dA-Za-z ]", "").replaceAll("\\s+", "+");
			   line = line.trim().replaceAll(" +", "");

					   if(line.length()>80)
							prot_id = line.substring(1, 81);
					   else
						   prot_id  = line.substring(1, line.length());
			
					   fastaString+=">"+prot_id+"\n";
			 
			 
		 if(lst_training_ids.contains(prot_id)){ 
			line = br.readLine();
		    while(line.startsWith(">")==false){//until > add sequence 
		    	 fastaString+=line;
				line = br.readLine();
				if(line==null)
					 break;
			   }
			fastaArray.put(prot_id,fastaString);
		 	} 
		 else{
				line = br.readLine();
			 while(line.startsWith(">")==false){//until > add sequence 
					line = br.readLine();
					if(line==null)
						 break;
				   }
		 }
		}	
		 
		for( i=0; i<lst_training_ids.size();i++){
			fastaWrite.add(fastaArray.get(lst_training_ids.get(i)));
		}
		for( i=0; i<fastaWrite.size();i++){
			if(fastaWrite.get(i)==null)
				continue;
			final_file.println(fastaWrite.get(i));
		}
		 br.close();
		 final_file.close();			 
	}
	
	
}
