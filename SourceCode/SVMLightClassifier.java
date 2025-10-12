import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

/**
 * In-JVM SVM classifier for SVMlight models.
 * Replaces external svm_classify process calls for better performance.
 * Optimized for speed with array-based support vectors and cached models.
 */
public class SVMLightClassifier {
    
    // Model cache to avoid reloading the same model (thread-safe)
    private static final Map<String, SVMLightClassifier> modelCache = new ConcurrentHashMap<>();
    
    private int kernelType;
    private int polyDegree;
    private double rbfGamma;
    private double coef0;
    private double threshold;
    private SupportVector[] supportVectors;  // Array for better performance
    private int maxFeatureIndex;
    
    private static class SupportVector {
        double alpha;  // alpha * y
        int[] indices;  // Feature indices (sorted)
        double[] values;  // Feature values
        
        SupportVector(double alpha, int size) {
            this.alpha = alpha;
            this.indices = new int[size];
            this.values = new double[size];
        }
    }
    
    /**
     * Load SVMlight model from file with caching
     */
    private void loadModel(String modelFile) throws IOException {
        List<SupportVector> svList = new ArrayList<>();
        maxFeatureIndex = 0;
        
        try (BufferedReader br = new BufferedReader(new FileReader(modelFile), 131072)) {
            // Line 1: SVM-light Version
            br.readLine();
            
            // Line 2: kernel type
            kernelType = Integer.parseInt(br.readLine().split("#")[0].trim());
            
            // Line 3: kernel parameter -d (degree for polynomial)
            polyDegree = Integer.parseInt(br.readLine().split("#")[0].trim());
            
            // Line 4: kernel parameter -g (gamma for RBF)
            rbfGamma = Double.parseDouble(br.readLine().split("#")[0].trim());
            
            // Line 5-8: other kernel parameters (s, r, u)
            br.readLine(); // -s
            coef0 = Double.parseDouble(br.readLine().split("#")[0].trim()); // -r
            br.readLine(); // -u
            
            // Line 9: highest feature index
            maxFeatureIndex = Integer.parseInt(br.readLine().split("#")[0].trim());
            
            // Line 10: number of training documents
            br.readLine();
            
            // Line 11: number of support vectors
            int numSVs = Integer.parseInt(br.readLine().split("#")[0].trim()) - 1;
            
            // Line 12: threshold b
            threshold = Double.parseDouble(br.readLine().split("#")[0].trim());
            
            // Pre-allocate array for support vectors
            svList = new ArrayList<>(numSVs);
            
            // Remaining lines: support vectors
            String line;
            while ((line = br.readLine()) != null) {
                int len = line.length();
                if (len == 0) continue;
                
                // Remove comments
                int commentIdx = line.indexOf('#');
                if (commentIdx >= 0) {
                    len = commentIdx;
                }
                
                // Fast parsing without split
                int pos = 0;
                
                // Skip leading whitespace
                while (pos < len && line.charAt(pos) == ' ') pos++;
                if (pos >= len) continue;
                
                // Parse alpha
                int spacePos = line.indexOf(' ', pos);
                if (spacePos == -1) continue;
                
                double alpha = Double.parseDouble(line.substring(pos, spacePos));
                pos = spacePos + 1;
                
                // Count features
                int numFeatures = 0;
                int tempPos = pos;
                while (tempPos < len) {
                    if (line.charAt(tempPos) == ':') numFeatures++;
                    tempPos++;
                }
                
                SupportVector sv = new SupportVector(alpha, numFeatures);
                
                // Parse features
                int featIdx = 0;
                while (pos < len && featIdx < numFeatures) {
                    // Skip whitespace
                    while (pos < len && line.charAt(pos) == ' ') pos++;
                    if (pos >= len) break;
                    
                    // Find colon
                    int colonPos = line.indexOf(':', pos);
                    if (colonPos == -1) break;
                    
                    // Parse index
                    sv.indices[featIdx] = parseInt(line, pos, colonPos);
                    pos = colonPos + 1;
                    
                    // Find next space or end
                    spacePos = line.indexOf(' ', pos);
                    if (spacePos == -1) spacePos = len;
                    
                    // Parse value
                    sv.values[featIdx] = Double.parseDouble(line.substring(pos, spacePos));
                    pos = spacePos;
                    featIdx++;
                }
                
                // Trim arrays if needed
                if (featIdx < numFeatures) {
                    sv.indices = Arrays.copyOf(sv.indices, featIdx);
                    sv.values = Arrays.copyOf(sv.values, featIdx);
                }
                
                svList.add(sv);
            }
        }
        
        // Convert to array for faster iteration
        supportVectors = svList.toArray(new SupportVector[0]);
    }
    
    /**
     * Fast integer parsing
     */
    private static int parseInt(String s, int start, int end) {
        int result = 0;
        boolean negative = false;
        int i = start;
        
        if (s.charAt(i) == '-') {
            negative = true;
            i++;
        }
        
        while (i < end) {
            result = result * 10 + (s.charAt(i) - '0');
            i++;
        }
        
        return negative ? -result : result;
    }
    
    /**
     * Classify a single instance (array-based for speed)
     */
    private double classify(int[] indices, double[] values) {
        double sum = 0.0;
        
        for (SupportVector sv : supportVectors) {
            double kernelValue = computeKernel(sv.indices, sv.values, indices, values);
            sum += sv.alpha * kernelValue;
        }
        
        return sum - threshold;
    }
    
    /**
     * Compute kernel function between two feature vectors (array-based)
     */
    private double computeKernel(int[] idx1, double[] val1, int[] idx2, double[] val2) {
        switch (kernelType) {
            case 0: // linear kernel
                return dotProduct(idx1, val1, idx2, val2);
                
            case 1: // polynomial kernel: (gamma * <x1,x2> + coef0)^degree
                double dot = dotProduct(idx1, val1, idx2, val2);
                return Math.pow(rbfGamma * dot + coef0, polyDegree);
                
            case 2: // RBF kernel: exp(-gamma * ||x1-x2||^2)
                return rbfKernel(idx1, val1, idx2, val2);
                
            case 3: // sigmoid kernel: tanh(gamma * <x1,x2> + coef0)
                dot = dotProduct(idx1, val1, idx2, val2);
                return Math.tanh(rbfGamma * dot + coef0);
                
            default:
                return dotProduct(idx1, val1, idx2, val2);
        }
    }
    
    /**
     * Compute dot product of two sparse vectors (optimized for sorted indices)
     */
    private double dotProduct(int[] idx1, double[] val1, int[] idx2, double[] val2) {
        double sum = 0.0;
        int i = 0, j = 0;
        
        // Both arrays are sorted, use merge-like algorithm
        while (i < idx1.length && j < idx2.length) {
            if (idx1[i] == idx2[j]) {
                sum += val1[i] * val2[j];
                i++;
                j++;
            } else if (idx1[i] < idx2[j]) {
                i++;
            } else {
                j++;
            }
        }
        
        return sum;
    }
    
    /**
     * Compute RBF kernel (optimized for sorted indices)
     */
    private double rbfKernel(int[] idx1, double[] val1, int[] idx2, double[] val2) {
        double squaredDist = 0.0;
        int i = 0, j = 0;
        
        // Handle elements in both vectors
        while (i < idx1.length && j < idx2.length) {
            if (idx1[i] == idx2[j]) {
                double diff = val1[i] - val2[j];
                squaredDist += diff * diff;
                i++;
                j++;
            } else if (idx1[i] < idx2[j]) {
                squaredDist += val1[i] * val1[i];
                i++;
            } else {
                squaredDist += val2[j] * val2[j];
                j++;
            }
        }
        
        // Handle remaining elements in first vector
        while (i < idx1.length) {
            squaredDist += val1[i] * val1[i];
            i++;
        }
        
        // Handle remaining elements in second vector
        while (j < idx2.length) {
            squaredDist += val2[j] * val2[j];
            j++;
        }
        
        return Math.exp(-rbfGamma * squaredDist);
    }
    
    /**
     * Classify all instances from a file and write predictions to output
     * Optimized version that doesn't store all instances in memory
     */
    private void classifyFile(String testFile, String predFile) throws IOException {
        // Pre-allocate StringBuilder for better performance
        StringBuilder lineBuffer = new StringBuilder(2048);
        
        try (BufferedReader br = new BufferedReader(new FileReader(testFile), 131072);
             BufferedWriter bw = new BufferedWriter(new FileWriter(predFile), 131072)) {
            
            String line;
            int[] indices = new int[maxFeatureIndex + 1];
            double[] values = new double[maxFeatureIndex + 1];
            
            while ((line = br.readLine()) != null) {
                int len = line.length();
                if (len == 0) continue;
                
                // Parse instance inline without split
                int featIdx = 0;
                int i = 0;
                
                // Skip label and first space
                while (i < len && line.charAt(i) != ' ') i++;
                i++; // skip space
                
                // Parse features
                while (i < len) {
                    // Read index
                    int colonPos = line.indexOf(':', i);
                    if (colonPos == -1) break;
                    
                    int index = 0;
                    while (i < colonPos) {
                        index = index * 10 + (line.charAt(i) - '0');
                        i++;
                    }
                    i++; // skip colon
                    
                    // Read value
                    int spacePos = line.indexOf(' ', i);
                    if (spacePos == -1) spacePos = len;
                    
                    double value = Double.parseDouble(line.substring(i, spacePos));
                    i = spacePos + 1;
                    
                    indices[featIdx] = index;
                    values[featIdx] = value;
                    featIdx++;
                }
                
                // Classify
                double prediction = classifyFast(indices, values, featIdx);
                
                // Write result
                bw.write(Double.toString(prediction));
                bw.write('\n');
            }
        }
    }
    
    /**
     * Fast double parsing for ASCII numeric strings
     */
    private static double parseDouble(String s, int start, int end) {
        double result = 0.0;
        boolean negative = false;
        int i = start;
        
        if (s.charAt(i) == '-') {
            negative = true;
            i++;
        }
        
        // Parse integer part
        while (i < end && s.charAt(i) != '.') {
            result = result * 10 + (s.charAt(i) - '0');
            i++;
        }
        
        // Parse fractional part
        if (i < end && s.charAt(i) == '.') {
            i++;
            double factor = 0.1;
            while (i < end) {
                result += (s.charAt(i) - '0') * factor;
                factor *= 0.1;
                i++;
            }
        }
        
        return negative ? -result : result;
    }
    
    /**
     * Fast classification without array allocation
     */
    private double classifyFast(int[] indices, double[] values, int numFeatures) {
        double sum = 0.0;
        
        // Unroll common cases
        if (kernelType == 0) {
            // Linear kernel - most common case
            for (SupportVector sv : supportVectors) {
                double dot = 0.0;
                int i = 0, j = 0;
                
                while (i < sv.indices.length && j < numFeatures) {
                    if (sv.indices[i] == indices[j]) {
                        dot += sv.values[i] * values[j];
                        i++;
                        j++;
                    } else if (sv.indices[i] < indices[j]) {
                        i++;
                    } else {
                        j++;
                    }
                }
                sum += sv.alpha * dot;
            }
        } else if (kernelType == 2) {
            // RBF kernel - second most common
            for (SupportVector sv : supportVectors) {
                double squaredDist = 0.0;
                int i = 0, j = 0;
                
                while (i < sv.indices.length && j < numFeatures) {
                    if (sv.indices[i] == indices[j]) {
                        double diff = sv.values[i] - values[j];
                        squaredDist += diff * diff;
                        i++;
                        j++;
                    } else if (sv.indices[i] < indices[j]) {
                        squaredDist += sv.values[i] * sv.values[i];
                        i++;
                    } else {
                        squaredDist += values[j] * values[j];
                        j++;
                    }
                }
                
                while (i < sv.indices.length) {
                    squaredDist += sv.values[i] * sv.values[i];
                    i++;
                }
                
                while (j < numFeatures) {
                    squaredDist += values[j] * values[j];
                    j++;
                }
                
                sum += sv.alpha * Math.exp(-rbfGamma * squaredDist);
            }
        } else {
            // Generic kernel computation
            for (SupportVector sv : supportVectors) {
                double kernelValue = computeKernelFast(sv.indices, sv.values, indices, values, numFeatures);
                sum += sv.alpha * kernelValue;
            }
        }
        
        return sum - threshold;
    }
    
    /**
     * Fast kernel computation
     */
    private double computeKernelFast(int[] idx1, double[] val1, int[] idx2, double[] val2, int len2) {
        switch (kernelType) {
            case 0: // linear kernel
                return dotProductFast(idx1, val1, idx2, val2, len2);
                
            case 1: // polynomial kernel
                double dot = dotProductFast(idx1, val1, idx2, val2, len2);
                return Math.pow(rbfGamma * dot + coef0, polyDegree);
                
            case 2: // RBF kernel
                return rbfKernelFast(idx1, val1, idx2, val2, len2);
                
            case 3: // sigmoid kernel
                dot = dotProductFast(idx1, val1, idx2, val2, len2);
                return Math.tanh(rbfGamma * dot + coef0);
                
            default:
                return dotProductFast(idx1, val1, idx2, val2, len2);
        }
    }
    
    private double dotProductFast(int[] idx1, double[] val1, int[] idx2, double[] val2, int len2) {
        double sum = 0.0;
        int i = 0, j = 0;
        
        while (i < idx1.length && j < len2) {
            if (idx1[i] == idx2[j]) {
                sum += val1[i] * val2[j];
                i++;
                j++;
            } else if (idx1[i] < idx2[j]) {
                i++;
            } else {
                j++;
            }
        }
        
        return sum;
    }
    
    private double rbfKernelFast(int[] idx1, double[] val1, int[] idx2, double[] val2, int len2) {
        double squaredDist = 0.0;
        int i = 0, j = 0;
        
        while (i < idx1.length && j < len2) {
            if (idx1[i] == idx2[j]) {
                double diff = val1[i] - val2[j];
                squaredDist += diff * diff;
                i++;
                j++;
            } else if (idx1[i] < idx2[j]) {
                squaredDist += val1[i] * val1[i];
                i++;
            } else {
                squaredDist += val2[j] * val2[j];
                j++;
            }
        }
        
        while (i < idx1.length) {
            squaredDist += val1[i] * val1[i];
            i++;
        }
        
        while (j < len2) {
            squaredDist += val2[j] * val2[j];
            j++;
        }
        
        return Math.exp(-rbfGamma * squaredDist);
    }
    
    /**
     * Static convenience method for classification with caching
     */
    public static void classify(String testFile, String modelFile, String predFile) throws IOException {
        // Check cache first (lock-free read for common case)
        SVMLightClassifier classifier = modelCache.get(modelFile);
        
        if (classifier == null) {
            // Only synchronize on cache miss
            classifier = new SVMLightClassifier();
            classifier.loadModel(modelFile);
            modelCache.putIfAbsent(modelFile, classifier);
            // Re-read in case another thread added it
            classifier = modelCache.get(modelFile);
        }
        
        classifier.classifyFile(testFile, predFile);
    }
    
    /**
     * Clear the model cache (useful for memory management)
     */
    public static void clearCache() {
        synchronized (modelCache) {
            modelCache.clear();
        }
    }
}
