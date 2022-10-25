import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * This program implements Brute Force Median Search, Branch and Bound Median Search and Greedy Motif Search algorithms
 * from the book Introduction to Bioinformatics Algorithms by Jones and Pevzner. It accepts file paths as command-line
 * arguments. Inputs files should have only 2 lines. The first line should have dna sequences of equal lengths separated
 * by commas and the second line should have an integer specifying the l (length of motif/l-mer). It outputs running
 * times for each of the three algorithms along with scores of their results over the corresponding DNA sequences. In
 * addition to processing inputs from files, numerous random DNA sequences are generated varying in one of n, t, or l at
 * a time. These sequences then go through the same process as inputs from files and the outputs and outputted to
 * console. For inputs from files, return values of methods are also outputted to console while for randomly generated
 * DNA sequences, they are not.
 */
class Code {
    public static void main(String[] args) {
        if (args.length > 0) {
            System.out.println("A. DNA Sequences from Files Paths Provided as Command-Line Arguments:");

            for (String arg : args) {
                Scanner inFile;
                try {
                    inFile = new Scanner(new File(arg));
                } catch (FileNotFoundException e) {
                    System.out.println("File " + arg + " was not found.\n");

                    continue;
                }

                String[] dna = inFile.nextLine().toUpperCase().split(",");
                int l = inFile.nextInt();

                inFile.close();

                runAlgorithmsAndPrint(dna, l, true);

                System.out.print("\n");
            }
        }

        System.out.println("B. DNA Sequences Generated at Random:");

        for (int n = 10; n <= 100 ; n += 10) {
            runAlgorithmsAndPrint(generateRandomDNASequences(5, n), 8, false);

            System.out.print("\n");
        }

        for (int t = 2; t <= 10 ; t += 2) {
            runAlgorithmsAndPrint(generateRandomDNASequences(t, 30), 8, false);

            System.out.print("\n");
        }

        for (int l = 5; l <= 10 ; l++) {
            runAlgorithmsAndPrint(generateRandomDNASequences(5, 30), l, false);

            System.out.print("\n");
        }
    }

    /**
     * Times the three algorithms (10 times each) and prints average times along with scores to console
     * @param dna receives DNA sequences
     * @param l length of motif/l-mer
     * @param printReturnVals if true, prints values returned by the three algorithms to console; if false, does not
     */
    private static void runAlgorithmsAndPrint(String[] dna, int l, boolean printReturnVals) {
        System.out.println("t = " + dna.length + "; n = " + dna[0].length() + "; l = " + l);

        long stTime, enTime;
        int repeatCount = 10;
        String stringRes = "";
        int[] arrayRes = {};

        stTime = System.nanoTime();
        for (int i = 0; i < repeatCount; i++) stringRes = bruteForceMedianSearch(dna, l);
        enTime = System.nanoTime();

        System.out.println("Brute Force Median Search; Average Time Taken: " + (enTime - stTime) / repeatCount +
                " nanoseconds; Score: " + (stringRes.length() * dna.length - totalDistance(stringRes, dna)) +
                (printReturnVals ? "; Return Value: " + stringRes : ""));

        stTime = System.nanoTime();
        for (int i = 0; i < repeatCount; i++) stringRes = branchAndBoundMedianSearch(dna, l);
        enTime = System.nanoTime();

        System.out.println("Branch and Bound Median Search; Average Time Taken: " + (enTime - stTime) / repeatCount +
                " nanoseconds; Score: " + (stringRes.length() * dna.length - totalDistance(stringRes, dna)) +
                (printReturnVals ? "; Return Value: " + stringRes : ""));

        stTime = System.nanoTime();
        for (int i = 0; i < repeatCount; i++) arrayRes = greedyMotifSearch(dna, l);
        enTime = System.nanoTime();

        System.out.println("Greedy Motif Search; Average Time Taken: " + (enTime - stTime) / repeatCount +
                " nanoseconds; Score: " + score(arrayRes, dna.length, dna, l) + (printReturnVals ? "; Return Value: " +
                Arrays.toString(arrayRes) : ""));
    }

    /**
     * Implementation of the algorithm mentioned on page 112 of Introduction to Bioinformatics Algorithms Book. Finds
     * the l-mer with the lowest total distance over 'dna'.
     * @param dna receives DNA sequences
     * @param l receives the length of l-mer
     * @return an l-mer of length l
     */
    private static String bruteForceMedianSearch(String[] dna, int l) {
        List<String> lMers = new ArrayList<>();
        addLMers(lMers, l, "");

        String bestWord = lMers.get(0);
        int bestDistance = totalDistance(lMers.get(0), dna);

        for (int i = 1; i < lMers.size(); i++) {
            if (totalDistance(lMers.get(i), dna) < bestDistance) {
                bestDistance = totalDistance(lMers.get(i), dna);
                bestWord = lMers.get(i);
            }
        }

        return bestWord;
    }

    /**
     * Implementation of the algorithm mentioned on page 114 of Introduction to Bioinformatics Algorithms Book. Finds
     * the l-mer with the lowest total distance over 'dna'.
     * @param dna receives DNA sequences
     * @param l receives the length of l-mer
     * @return an l-mer of length l
     */
    private static String branchAndBoundMedianSearch(String[] dna, int l) {
        int[] s = new int[l];
        s[0] = 1;

        String bestWord = "";
        int bestDistance = Integer.MAX_VALUE;

        int i = 1;

        while (i > 0) {
            if (i < l) {
                String prefix = intArrayToNucleotideString(s, i);
                int optimisticDistance = totalDistance(prefix, dna);

                if (optimisticDistance > bestDistance) {
                    i = bypass(s, i);
                } else {
                    i = nextVertex(s, i, l);
                }
            } else {
                String word = intArrayToNucleotideString(s, l);

                if (totalDistance(word, dna) < bestDistance) {
                    bestDistance = totalDistance(word, dna);
                    bestWord = word;
                }

                i = nextVertex(s, i, l);
            }
        }

        return bestWord;
    }

    /**
     * Implementation of the algorithm mentioned on page 136 of Introduction to Bioinformatics Algorithms Book. Attempts
     * to find motifs in sequences in 'dna' that would lead to a score close to the highest possible score.
     * @param dna receives DNA sequences
     * @param l receives the length of motif
     * @return indices of start positions for each DNA sequence in 'dna'
     */
    private static int[] greedyMotifSearch(String[] dna, int l) {
        int[] bestMotif = new int[dna.length];
        Arrays.fill(bestMotif, 0);

        int[] s = new int[dna.length];
        Arrays.fill(s, 0);

        for (int i = 0; i < dna[0].length() - l + 1; i++) {
            s[0] = i;

            for (int j = 0; j < dna[0].length() - l + 1; j++) {
                s[1] = j;

                if (score(s, 2, dna, l) > score(bestMotif, 2, dna, l)) {
                    bestMotif[0] = s[0];
                    bestMotif[1] = s[1];
                }
            }
        }

        s[0] = bestMotif[0];
        s[1] = bestMotif[1];

        for (int i = 2; i < dna.length; i++) {
            for (int j = 0; j < dna[i].length() - l + 1; j++) {
                s[i] = j;

                if (score(s, i + 1, dna, l) > score(bestMotif, i + 1, dna, l)) bestMotif[i] = s[i];
            }

            s[i] = bestMotif[i];
        }

        return bestMotif;
    }

    /**
     * Creates and adds all combinations of A, T, G and C of length 'l' to 'words'
     * @param words receives the list to which l-mers are added
     * @param l receives length of each l-mer
     * @param word an l-mer that recursively changes on each iteration
     */
    private static void addLMers(List<String> words, int l, String word) {
        if (word.length() == l) {
            words.add(word);
        } else {
            addLMers(words, l, word + "A");
            addLMers(words, l, word + "T");
            addLMers(words, l, word + "G");
            addLMers(words, l, word + "C");
        }
    }

    /**
     * Calculates the total distance of the 'word' over the 'dna' sequences
     * @param word receives a combination og nucleotides A, T, G and C
     * @param dna receives DNA sequences
     * @return the calculated total distance
     */
    private static int totalDistance(String word, String[] dna) {
        int[] dnaLowestDistances = new int[dna.length];
        Arrays.fill(dnaLowestDistances, -1);

        for (int i = 0; i < dna.length; i++) {
            for (int j = 0; j < dna[i].length() - word.length() + 1; j++) {
                int distance = 0;
                for (int k = 0; k < word.length(); k++) if (word.charAt(k) != dna[i].charAt(j + k)) distance += 1;

                if (dnaLowestDistances[i] == -1 || distance < dnaLowestDistances[i]) dnaLowestDistances[i] = distance;
            }
        }

        int totalDistance = 0;
        for (int dnaLowestDistance : dnaLowestDistances) totalDistance += dnaLowestDistance;

        return totalDistance;
    }

    /**
     * Converts an integer array to a DNA sequence string containing nucleotides
     * @param s receives the representation a DNA sequence in digits; 1 = A; 2 = T; 3 = G; 4 = C
     * @param i receives length of the nucleotide string
     * @return the corresponding nucleotide string
     */
    private static String intArrayToNucleotideString(int[] s, int i) {
        StringBuilder nucleotideString = new StringBuilder();

        for (int j = 0; j < i; j++) {
            if (s[j] == 1) {
                nucleotideString.append("A");
            } else if (s[j] == 2) {
                nucleotideString.append("T");
            } else if (s[j] == 3) {
                nucleotideString.append("G");
            } else {
                nucleotideString.append("C");
            }
        }

        return nucleotideString.toString();
    }

    /**
     * Implementation of the algorithm mentioned on page 107 of Introduction to Bioinformatics Algorithms Book. Updates
     * 's' to hold the next prefix/vertex and returns its length.
     * @param s receives the representation a DNA sequence in digits; 1 = A; 2 = T; 3 = G; 4 = C
     * @param i receives length of prefix/vertex
     * @param l receives maximum length
     * @return length of next prefix/vertex; 0 if end reached
     */
    private static int nextVertex(int[] s, int i, int l) {
        if (i < l) {
            s[i] = 1;
            return i + 1;
        } else {
            for (int j = l; j >= 1; j--) {
                if (s[j - 1] < 4) {
                    s[j - 1]++;
                    return j;
                }
            }
        }

        return 0;
    }

    /**
     * Implementation of the algorithm mentioned on page 108 of Introduction to Bioinformatics Algorithms Book. Updates
     * 's' (skipping children) to hold the next sibling prefix/vertex and returns its length.
     * @param s receives the representation a DNA sequence in digits; 1 = A; 2 = T; 3 = G; 4 = C
     * @param i receives length of prefix/vertex
     * @return length of next prefix/vertex; 0 if end reached
     */
    private static int bypass(int[] s, int i) {
        for (int j = i; j >= 1; j--) {
            if (s[j - 1] < 4) {
                s[j - 1]++;
                return j;
            }
        }

        return 0;
    }

    /**
     * Calculates score of first 'numOfSeq' sequences in 'dna'
     * @param s receives indices of start positions for each DNA sequence in 'dna'
     * @param numOfSeq receives the number of sequences to score
     * @param dna receives DNA sequences
     * @param l receives the length of motif
     * @return the calculated score
     */
    private static int score(int[] s, int numOfSeq, String[] dna, int l) {
        int score = 0;

        for (int i = 0; i < l; i++) {
            Integer[] counts = {0, 0, 0, 0};

            for (int j = 0; j < numOfSeq; j++) {
                switch (dna[j].charAt(s[j] + i)) {
                    case 'A':
                        counts[0]++;
                        break;
                    case 'T':
                        counts[1]++;
                        break;
                    case 'G':
                        counts[2]++;
                        break;
                    case 'C':
                        counts[3]++;
                }
            }

            score += Collections.max(Arrays.asList(counts));
        }

        return score;
    }

    /**
     * Generates random DNA sequences of nucleotides conforming to received specifications
     * @param numOfSeq receives the number of sequences to generate
     * @param lenOfEachSeq receives the length of each sequence
     * @return the generated sequences
     */
    private static String[] generateRandomDNASequences(int numOfSeq, int lenOfEachSeq) {
        String[] dna = new String[numOfSeq];
        String[] nucleotides = {"A", "T", "G", "C"};

        for (int i = 0; i < numOfSeq; i++) {
            StringBuilder sequence = new StringBuilder();

            for (int j = 0; j < lenOfEachSeq; j++) sequence.append(nucleotides[new Random().nextInt(4)]);

            dna[i] = sequence.toString();
        }

        return dna;
    }
}