# Bug Fix: ArrayIndexOutOfBoundsException for Short Proteins

## Problem Description

The program crashes with an `ArrayIndexOutOfBoundsException` when processing protein sequences that are too short (less than 5 amino acids):

```
Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: Array index out of range: 1545
        at java.base/java.util.Vector.get(Vector.java:750)
        at runEC.processMainClassPredictions(runEC.java:241)
        at runEC.predictions(runEC.java:81)
        at ECPred.main(ECPred.java:130)
```

## Root Cause

The issue occurs because of a mismatch between the number of proteins in the input and the number of predictions generated:

### 1. **seq2vectPSSMtest.java** (SPMAP method)
- **Line 88-89**: Short proteins (< 5 amino acids) are skipped with `continue`
- **Result**: No feature vector is generated for these proteins
- **Impact**: Prediction files have fewer lines than proteins

### 2. **predictBatchBLAST.java** (BLAST method)
- **Line 143**: When BLAST has no hits, only ONE 0.0 value is added instead of one per protein
- **Result**: Prediction file has only 1 line when it should have N lines
- **Impact**: Severe mismatch in prediction count

### 3. **fasta2Pepstats_noscale.java** (PEPSTATS method)
- **Line 92-98**: Proteins with "None" in pepstats output are skipped
- **Result**: No feature vector is generated for these proteins
- **Impact**: Prediction files have fewer lines than proteins

### 4. **runEC.java** - Where the crash occurs
- **Line 241**: `allPreds.get(j).get(i)` assumes all proteins have predictions
- Iterates `i` from 0 to `idlist.size()-1` (all proteins)
- Tries to access prediction by index, but prediction vectors are shorter
- **Result**: ArrayIndexOutOfBoundsException when accessing skipped proteins

## The Fix

All three prediction methods have been fixed to ensure they generate predictions for ALL proteins, even if some are too short or cannot be processed:

### Fix 1: seq2vectPSSMtest.java (Lines 74-107)

**Before:**
```java
for (String targetId : targetList) {
    String sequence = fastaDict.get(targetId);
    if (sequence == null || sequence.length() < subseqlen) {
        continue;  // SKIPS the protein
    }
    
    double[] bestScores = computeBestScores(sequence, profile, subseqlen);
    if (bestScores == null) {
        continue;  // SKIPS the protein
    }
    
    applySignificanceThreshold(bestScores, signifThreshold, subseqlen);
    vectors.add(bestScores);
}
```

**After:**
```java
int clusterCount = profile.weights.length;

for (String targetId : targetList) {
    String sequence = fastaDict.get(targetId);
    if (sequence == null || sequence.length() < subseqlen) {
        // Add a zero vector for short/missing sequences
        double[] zeroVector = new double[clusterCount];
        Arrays.fill(zeroVector, 0.0);
        vectors.add(zeroVector);
        continue;
    }
    
    double[] bestScores = computeBestScores(sequence, profile, subseqlen);
    if (bestScores == null) {
        // Add a zero vector for sequences that couldn't be processed
        double[] zeroVector = new double[clusterCount];
        Arrays.fill(zeroVector, 0.0);
        vectors.add(zeroVector);
        continue;
    }
    
    applySignificanceThreshold(bestScores, signifThreshold, subseqlen);
    vectors.add(bestScores);
}
```

### Fix 2: predictBatchBLAST.java (Lines 137-156)

**Before:**
```java
List<String> blastLines = Files.readAllLines(Paths.get(cmdArray[8], new String[0]));
if (blastLines.size() == 0) {
    preds.add(Double.valueOf(0.0D));  // Only ONE prediction!
} else {
    simHashHash = Blast.parseTabBlast(cmdArray[8]);
    preds = new Vector<>();
    for (int m = 0; m < test_ids.size(); m++) {
        double pred = Blast.blastknn(simHashHash.get(test_ids.get(m)), pos, neg, k);
        preds.add(Double.valueOf(pred));
    }
}
```

**After:**
```java
List<String> blastLines = Files.readAllLines(Paths.get(cmdArray[8], new String[0]));
if (blastLines.size() == 0) {
    // Add 0.0 prediction for EACH test protein when BLAST has no hits
    for (int m = 0; m < test_ids.size(); m++) {
        preds.add(Double.valueOf(0.0D));
    }
} else {
    simHashHash = Blast.parseTabBlast(cmdArray[8]);
    preds = new Vector<>();
    for (int m = 0; m < test_ids.size(); m++) {
        double pred = Blast.blastknn(simHashHash.get(test_ids.get(m)), pos, neg, k);
        preds.add(Double.valueOf(pred));
    }
}
```

### Fix 3: fasta2Pepstats_noscale.java (Lines 74-108)

**Before:**
```java
while ((line = br.readLine()) != null) {
    count++;
    single.add(line);
    if (line.contains("None")) {
        while (count % 49 != 0) {
            line = br.readLine();
            count++;
        }
        count = 0;
        single = new Vector<>();  // SKIPS this protein
    } else if (count % 48 == 0) {
        vects.add(vecCount, parse_single(single));
        vecCount++;
        single = new Vector<>();
        count = 0;
    }
}
```

**After:**
```java
while ((line = br.readLine()) != null) {
    count++;
    single.add(line);
    if (line.contains("None")) {
        while (count % 49 != 0) {
            line = br.readLine();
            count++;
        }
        // Add a zero vector for proteins with "None"
        Vector<String> zeroVector = new Vector<>();
        for (int i = 0; i < 40; i++) { // pepstats has 40 features
            zeroVector.add("0.0");
        }
        vects.add(vecCount, zeroVector);
        vecCount++;
        count = 0;
        single = new Vector<>();
    } else if (count % 48 == 0) {
        vects.add(vecCount, parse_single(single));
        vecCount++;
        single = new Vector<>();
        count = 0;
    }
}
```

## How the Fix Works

Instead of skipping proteins that are too short or cannot be processed:

1. **Generate a zero/default feature vector** for problematic proteins
2. **Maintain the 1-to-1 correspondence** between input proteins and predictions
3. **Ensure prediction files have the same number of lines** as proteins in the input

This way:
- Short proteins get a prediction score of 0.0 (classified as non-enzyme)
- The indexing in `runEC.processMainClassPredictions()` works correctly
- No ArrayIndexOutOfBoundsException occurs

## Testing

After applying these fixes, the program should:
- Successfully process FASTA files containing short proteins
- Generate predictions for all proteins in the input
- No longer crash with ArrayIndexOutOfBoundsException

## Impact

- **Compatibility**: No changes to existing functionality for normal-length proteins
- **Behavior**: Short proteins will now be classified as non-enzymes (score 0.0) instead of causing crashes
- **Robustness**: Program can now handle mixed-length protein datasets
