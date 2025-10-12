import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.concurrent.ConcurrentHashMap;

public class seq2vectPSSMtest
{
  private static final int MAX_PROTEIN_ID_LENGTH = 80;
  private static final int MIN_SUBSEQUENCE_COUNT = 6;
  private static final int[] AA_INDEX = new int[26];
  private static final int AA_COUNT = 20;
  private static final String AA_ORDER = "ARNDCQEGHILKMFPSTWYV";
  private static final ConcurrentHashMap<String, Map<String, String>> FASTA_CACHE = new ConcurrentHashMap<>();
  private static final ConcurrentHashMap<String, PSSMProfile> PSSM_CACHE = new ConcurrentHashMap<>();
  
  static {
    Arrays.fill(AA_INDEX, -1);
    for (int i = 0; i < AA_ORDER.length(); i++) {
      AA_INDEX[AA_ORDER.charAt(i) - 'A'] = i;
    }
  }
  
  public static HashMap<String, String> getFastaOrg(String file)
    throws IOException
  {
    HashMap<String, String> fasta_dict = new HashMap<>();
    String prot_id = "";
    String prot_seq = "";
    
    BufferedReader br = new BufferedReader(new FileReader(file));
    String line = br.readLine();
    while (line != null)
    {
      if (line.startsWith(">"))
      {
        prot_seq = "";
        StringTokenizer st1 = new StringTokenizer(line, "|");
        line = line.replaceAll("/", "");
        line = line.trim().replaceAll(" +", "");
        if (line.length() > MAX_PROTEIN_ID_LENGTH) {
          prot_id = line.substring(1, MAX_PROTEIN_ID_LENGTH + 1);
        } else {
          prot_id = line.substring(1, line.length());
        }
        line = br.readLine();
        while (!line.startsWith(">"))
        {
          prot_seq = prot_seq + line;
          line = br.readLine();
          if (line == null) {
            break;
          }
        }
      }
      fasta_dict.put(prot_id, prot_seq);
    }
    br.close();
    return fasta_dict;
  }
  
  public static float calculateVectors(int signifThreshold, int subseqlen, String ECNumber, List<String> targetList, String filename, long time, String ROOTPATH, String tempDir)
    throws IOException
  {
    Map<String, String> fastaDict = getCachedFasta(filename);
    PSSMProfile profile = getCachedPSSMProfile(ROOTPATH, ECNumber, subseqlen);
    if (profile == null) {
      return 0.0F;
    }

    List<double[]> vectors = new ArrayList<>();
    for (String targetId : targetList)
    {
      String sequence = fastaDict.get(targetId);
      if (sequence == null || sequence.length() < subseqlen) {
        continue;
      }

      double[] bestScores = computeBestScores(sequence, profile, subseqlen);
      if (bestScores == null) {
        continue;
      }

      applySignificanceThreshold(bestScores, signifThreshold, subseqlen);
      vectors.add(bestScores);
    }

    writeVectors(ECNumber, vectors, time, tempDir);
    return 0.0F;
  }
  
  private static Map<String, String> getCachedFasta(String filename) throws IOException {
    Map<String, String> cached = FASTA_CACHE.get(filename);
    if (cached != null) {
      return cached;
    }
    Map<String, String> loaded = Collections.unmodifiableMap(getFastaOrg(filename));
    Map<String, String> existing = FASTA_CACHE.putIfAbsent(filename, loaded);
    return existing != null ? existing : loaded;
  }

  private static PSSMProfile getCachedPSSMProfile(String ROOTPATH, String ECNumber, int subseqlen) throws IOException {
    PSSMProfile cached = PSSM_CACHE.get(ECNumber);
    if (cached != null) {
      return cached;
    }
    PSSMProfile profile = loadPSSMProfile(ROOTPATH, ECNumber, subseqlen);
    if (profile == null) {
      return null;
    }
    PSSMProfile existing = PSSM_CACHE.putIfAbsent(ECNumber, profile);
    return existing != null ? existing : profile;
  }

  private static PSSMProfile loadPSSMProfile(String ROOTPATH, String ECNumber, int subseqlen)
    throws IOException
  {
    Path profilePath = Paths.get(ROOTPATH + File.separator + ECNumber + "/spmap/profile.txt", new String[0]);
    if (!Files.exists(profilePath)) {
      return null;
    }

    List<String> lines = Files.readAllLines(profilePath);
    if (lines.isEmpty()) {
      return null;
    }

    int numberOfCluster = lines.size() / subseqlen;
    double[][][] weights = new double[numberOfCluster][subseqlen][AA_COUNT];
    int count = 0;

    for (int cluster = 0; cluster < numberOfCluster; cluster++)
    {
      for (int position = 0; position < subseqlen; position++)
      {
        String line = lines.get(count++);
        StringTokenizer st1 = new StringTokenizer(line, "\t");
        String list = null;
        while (st1.hasMoreElements()) {
          list = st1.nextToken();
        }
        if (list == null) {
          continue;
        }
        String[] aaTokens = list.split("\\s*,\\s*");
        for (int idx = 0; idx + 1 < aaTokens.length; idx += 2) {
          int aaIndex = toAminoAcidIndex(aaTokens[idx]);
          if (aaIndex >= 0) {
            weights[cluster][position][aaIndex] = parseDoubleSafe(aaTokens[idx + 1]);
          }
        }
      }
    }
    return new PSSMProfile(weights);
  }

  private static double parseDoubleSafe(String value) {
    try {
      return Double.parseDouble(value);
    } catch (NumberFormatException ex) {
      return 0.0D;
    }
  }

  private static int toAminoAcidIndex(String token) {
    if (token == null || token.isEmpty()) {
      return -1;
    }
    char c = Character.toUpperCase(token.charAt(0));
    if (c < 'A' || c > 'Z') {
      return -1;
    }
    return AA_INDEX[c - 'A'];
  }

  private static double[] computeBestScores(String sequence, PSSMProfile profile, int subseqlen) {
    int clusterCount = profile.weights.length;
    if (clusterCount == 0 || sequence.length() < subseqlen || sequence.length() < MIN_SUBSEQUENCE_COUNT) {
      return null;
    }

    double[] bestScores = new double[clusterCount];
    Arrays.fill(bestScores, Double.NEGATIVE_INFINITY);
    char[] chars = sequence.toCharArray();
    int[] aaIdxBuffer = new int[subseqlen];

    for (int start = 0; start <= chars.length - subseqlen; start++)
    {
      boolean valid = true;
      for (int offset = 0; offset < subseqlen; offset++)
      {
        char residue = Character.toUpperCase(chars[start + offset]);
        if (residue < 'A' || residue > 'Z') {
          valid = false;
          break;
        }
        int aaIndex = AA_INDEX[residue - 'A'];
        if (aaIndex < 0) {
          valid = false;
          break;
        }
        aaIdxBuffer[offset] = aaIndex;
      }

      if (!valid) {
        continue;
      }

      for (int cluster = 0; cluster < clusterCount; cluster++)
      {
        double sum = 0.0D;
        double[][] clusterWeights = profile.weights[cluster];
        for (int position = 0; position < subseqlen; position++)
        {
          sum += clusterWeights[position][aaIdxBuffer[position]];
        }
        if (sum > bestScores[cluster]) {
          bestScores[cluster] = sum;
        }
      }
    }

    boolean hasValid = false;
    for (double score : bestScores) {
      if (score != Double.NEGATIVE_INFINITY) {
        hasValid = true;
        break;
      }
    }
    if (!hasValid) {
      return null;
    }
    return bestScores;
  }

  private static void applySignificanceThreshold(double[] scores, int signifThreshold, int subseqlen) {
    for (int i = 0; i < scores.length; i++) {
      double value = scores[i];
      if (value == Double.NEGATIVE_INFINITY || value < signifThreshold) {
        scores[i] = 0.0D;
      } else {
        scores[i] = Math.exp(value / subseqlen);
      }
    }
  }

  public static void writeVectors(String ECNumber, List<double[]> result, long time, String tempDir)
    throws IOException
  {
    Path outputPath = Paths.get(tempDir, "testResult", String.valueOf(time), ECNumber, "spmap", "test.vec");
    Files.createDirectories(outputPath.getParent());
    try (BufferedWriter writer = Files.newBufferedWriter(outputPath, StandardCharsets.UTF_8)) {
      for (double[] vector : result)
      {
        StringBuilder sb = new StringBuilder("1");
        for (int key = 0; key < vector.length; key++) {
          sb.append(' ').append(key + 1).append(':').append(vector[key]);
        }
        writer.write(sb.toString());
        writer.newLine();
      }
    }
  }

  private static class PSSMProfile {
    private final double[][][] weights;

    private PSSMProfile(double[][][] weights) {
      this.weights = weights;
    }
  }
  
  public static HashMap<String, Integer> readBLOSUM62Matrix()
    throws IOException
  {
    String aa_letters = "A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V";
    String[] lst_aa_letters = aa_letters.split(",");
    HashMap<String, Integer> blosum_dict = new HashMap<>();
    BufferedReader br = new BufferedReader(new FileReader("blo62.csv"));
    String[] blo62_mat = null;
    
    int i = 0;
    String line;
    while ((line = br.readLine()) != null)
    {

      blo62_mat = line.split(",");
      for (int j = 0; j < blo62_mat.length; j++) {
        blosum_dict.put(lst_aa_letters[i] + "," + lst_aa_letters[j], Integer.valueOf(Integer.parseInt(blo62_mat[j])));
      }
      i++;
    }
    return blosum_dict;
  }
}