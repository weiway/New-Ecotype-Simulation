/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013  Jason M. Wood, Montana State University
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package ecosim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;


/**
 *  Run the demarcation confidence interval program.
 *
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class DemarcationConfidenceInterval {

    /**
     *  Run the demarcation confidence interval program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object for the entire dataset.
     *  @param samplePhylogeny The Phylogeny object of the sample. 
     *  @param binning The Binning object.
     *  @param hillclimb The Hillclimb object.
     */
    public DemarcationConfidenceInterval(MasterVariables masterVariables,
        Phylogeny phylogeny, Phylogeny samplePhylogeny, Binning binning,
        ParameterSet hillclimb) {
        this(
            masterVariables, phylogeny, samplePhylogeny, binning, hillclimb, ""
        );
    }

    /**
     *  Run the demarcation confidence interval program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object for the entire dataset.
     *  @param samplePhylogeny The Phylogeny object of the sample. 
     *  @param binning The Binning object.
     *  @param hillclimb The Hillclimb object.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public DemarcationConfidenceInterval(MasterVariables masterVariables,
        Phylogeny phylogeny, Phylogeny samplePhylogeny, Binning binning,
        ParameterSet hillclimb, String suffix) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;
        this.samplePhylogeny = samplePhylogeny;
        this.binning = binning;
        this.hillclimb = hillclimb;
        String workingDirectory = masterVariables.getWorkingDirectory();
        inputFileName = workingDirectory + "ES2demarcationIn" + suffix + ".dat";
        outputFileName = workingDirectory + "ES2demarcationOut" + suffix + ".dat";
        inputFileName2 = workingDirectory + "ES1demarcationIn" + suffix + ".dat";
        outputFileName2 = workingDirectory + "ES1demarcationOut" + suffix + ".dat";
        hasRun = false;
    }

    /**
     *  Run the demarcation confidence interval program.
     */
    public void run() {
        Execs execs = masterVariables.getExecs();
        int hashCode = this.hashCode();
        File inputFile = new File(inputFileName);
        File outputFile = new File(outputFileName);
        // Write the input values for the program to the demarcationIn.dat file.
		System.out.println("writing input file for ES2 demarcationsCI.app");
		System.out.println("Input File:");
		System.out.println(inputFileName);
        writeInputFile(inputFile);
        // Run the demarcationCI program.
        execs.runDemarcationCI(inputFile, outputFile);
        // Get the output provided by the demarcationCI program.
        readOutputFile(outputFile);
        // Set the flag stating that the confidence interval program has run.
        if (result > 0) {
            hasRun = true;
        }
    }
	
    public int[] run2() {
        Execs execs = masterVariables.getExecs();
        int hashCode = this.hashCode();
        File inputFile = new File(inputFileName2);
        File outputFile = new File(outputFileName2);
        // Write the input values for the program to the demarcationIn.dat file.
		System.out.println("writing input file for ES1 demarcationsCI.app");
		System.out.println("Input File:");
		System.out.println(inputFileName2);
        writeInputFile2(inputFile);
		System.out.println("Finished writing input file, now runing demarcationsCI.app");
        // Run the demarcationCI program.
		final long startTime = System.currentTimeMillis();
        execs.runDemarcationsCI(inputFile, outputFile);
		final long endTime = System.currentTimeMillis();
		System.out.println("Finished running demarcationsCI.app");
		System.out.println("Total execution time: " + (endTime - startTime)/1000 + "secs" );
		System.out.println("    ");
        // Get the output provided by the demarcationCI program.
        int [] interval = getDemarcOutput();
		// Set the flag stating that the confidence interval program has run.
        if (result > 0) {
            hasRun = true;
        }
		return interval;
    }
	

    /**
     *  Returns true if the demarcation confidence interval has been run, false
     *  otherwise.
     *
     *  @return True if the demarcation confidence interval has been run, false
     *  otherwise.
     */
    public boolean hasRun() {
        return hasRun;
    }

    /**
     *  Changes the value of hasRun.
     *
     *  @param hasRun The new value of hasRun.
     */
    public void setHasRun(boolean hasRun) {
        this.hasRun = hasRun;
    }

    /**
     *  Returns the result of the demarcation confidence interval program.
     *
     *  @return The result.
     */
    public int getResult() {
        return result;
    }

    /**
     *  Returns the likelihood of the result of the demarcation confidence
     *  interval program.
     *
     *  @return The likelihood.
     */
    public double getLikelihood() {
        return likelihood;
    }

    /**
     *  Returns the demarcation confidence interval as a String.
     *
     *  @return the demarcation confidence interval.
     */
    public String toString() {
        int npop = hillclimb.getNpop();
        return String.format("%d (%d)", npop, result);
    }
    /**
     *  Private method to write the input file for the demarcation confidence
     *  interval program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeInputFile(File inputFile) {
        ArrayList<BinLevel> bins = binning.getBinLevels();
        ParameterSet hillclimbResult = hillclimb;
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(inputFile));
            writer.write(String.format("%-20d numcrit\n", bins.size()));
			System.out.println("===============================================");
			System.out.print(String.format("%-20d numcrit\n", bins.size()));
            // Output just the number of bins at each crit level.
            for (int i = 0; i < bins.size(); i ++) {
                writer.write(String.format("%-20d\n", bins.get(i).getLevel()));
				System.out.print(String.format("%-20d\n", bins.get(i).getLevel()));
            }
            // Output the crit levels and the number of bins.
            for (int j = 0; j < bins.size(); j ++) {
                writer.write(String.format(
                    "%-20.6f %-20d\n",
                    bins.get(j).getCrit(),
                    bins.get(j).getLevel()
                ));
				System.out.print(String.format(
                    "%-20.6f %-20d\n",
                    bins.get(j).getCrit(),
                    bins.get(j).getLevel()
                ));
            }
            // Write the omega value.
            writer.write(
                String.format("%-20.5f omega\n", hillclimbResult.getOmega())
            );
			System.out.print(String.format("%-20.5f omega\n", hillclimbResult.getOmega()));
            // Write the sigma value.
            writer.write(
                String.format("%-20.5f sigma\n", hillclimbResult.getSigma())
            );
			System.out.print(String.format("%-20.5f sigma\n", hillclimbResult.getSigma()));
            // Estimate the value of npop by multiplying the npop value found
            // in hillclimbing by the ratio of the sample size to the total
            // number of environmental sequences.
            int npop = hillclimbResult.getNpop() * 
                samplePhylogeny.getNu() / phylogeny.getNu();
            if (npop < 1) {
                npop = 1;
            }
            // Write the npop value.
            writer.write(
                String.format("%-20d npop\n", npop)
            );
			System.out.print(String.format("%-20d npop\n", npop));
            // Write the step value.
            writer.write(
                String.format("%-20d step\n", masterVariables.getStep())
            );
			System.out.print(String.format("%-20d step\n", masterVariables.getStep()));
            // Write the nu value.
            writer.write(
                String.format("%-20d nu\n", samplePhylogeny.getNu())
            );
			System.out.print(String.format("%-20d nu\n", samplePhylogeny.getNu()));
            // Write the nrep value.
            writer.write(
                String.format("%-20d nrep\n", masterVariables.getNrep())
            );
			System.out.print(String.format("%-20d nrep\n", masterVariables.getNrep()));
            // Create the random number seed; an odd integer less than nine
            // digits long.
            long iii = (long)(100000000 * Math.random());
            if (iii % 2 == 0) {
                iii ++;
            }
            // Write the random number seed.
            writer.write(
                String.format("%-20d iii (random number seed)\n", iii)
            );
			System.out.print(String.format("%-20d iii (random number seed)\n", iii));
            // Write the length of the sequences.
            writer.write(
                String.format(
                    "%-20d lengthseq (after deleting gaps, etc.)\n",
                    samplePhylogeny.length()
                )
            );
			System.out.print(String.format("%-20d lengthseq (after deleting gaps, etc.)\n", samplePhylogeny.length()));
            // Write the whichavg value.
            int whichavg = masterVariables.getCriterion();
            writer.write(String.format("%-20d whichavg\n", whichavg));
			System.out.print(String.format("%-20d whichavg\n", whichavg));
            // Write the likelihoodsolution value.
            writer.write(
                String.format(
                    "%-20.5f likelihoodsolution\n",
                    hillclimbResult.getValue()
                )
            );
			System.out.print(String.format(
                    "%-20.5f likelihoodsolution\n",
                    hillclimbResult.getValue()
                )
			);
			System.out.println("===========================================");
        }
        catch (IOException e) {
            System.out.println("Error writing the input file for the " +
                               "confidence interval program.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the input file for the " +
                                       "confidence interval program.");
                }
            }
        }
    }

    private void writeInputFile2(File inputFile) {
        ArrayList<BinLevel> bins = binning.getBinLevels();
        ParameterSet hillclimbResult = hillclimb;
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(inputFile));
			System.out.println("------------------------------------------------------");
            writer.write(String.format("%-20d numcrit\n", bins.size()));
			System.out.print(String.format("%-20d numcrit\n", bins.size()));
            // Output just the number of bins at each crit level.
            for (int i = 0; i < bins.size(); i ++) {
                writer.write(String.format("%-20d\n", bins.get(i).getLevel()));
				System.out.print(String.format("%-20d\n", bins.get(i).getLevel()));
            }
            // Output the crit levels and the number of bins.
            for (int j = 0; j < bins.size(); j ++) {
                writer.write(String.format(
                    "%-20.6f %-20d\n",
                    bins.get(j).getCrit(),
                    bins.get(j).getLevel()
                ));
				System.out.print(String.format(
                    "%-20.6f %-20d\n",
                    bins.get(j).getCrit(),
                    bins.get(j).getLevel()
                ));
            }
            // Write the omega value.
            writer.write(String.valueOf(hillclimbResult.getOmega()) + ",10000.0         omega\n");
			System.out.print(String.valueOf(hillclimbResult.getOmega()) + ",10000.0         omega\n");
			// Write the sigma value.
            writer.write(String.valueOf(hillclimbResult.getSigma()) + ",10000.0         sigma\n");
			System.out.print(String.valueOf(hillclimbResult.getSigma()) + ",10000.0         sigma\n");
			// Estimate the value of npop by multiplying the npop value found
            // in hillclimbing by the ratio of the sample size to the total
            // number of environmental sequences.
            int npop = hillclimbResult.getNpop() * 
                samplePhylogeny.getNu() / phylogeny.getNu();
            if (npop < 1) {
                npop = 1;
            }
            // Write the npop value.
            writer.write("1,3                npop\n");
			System.out.print("1,3                npop\n");
			// Write drift value.
			writer.write("1.0E25,1.0E26        drift\n");
			System.out.print("1.0E25,1.0E26        drift\n");
            // Write the xnumincs value.
            writer.write("0,0,1,0              xnumincs\n");
			System.out.print("0,0,1,0              xnumincs\n");
            // Write the nu value.
            writer.write(
                String.format("%-20d nu\n", samplePhylogeny.getNu())
            );
			System.out.print(String.format("%-20d nu\n", samplePhylogeny.getNu()));
            // Write the nrep value.
            writer.write(
                String.format("500                  nrep\n")
            );
			System.out.print(String.format("500                  nrep\n"));
            // Create the random number seed; an odd integer less than nine
            // digits long.
            long iii = (long)(100000000 * Math.random());
            if (iii % 2 == 0) {
                iii ++;
            }
            // Write the random number seed.
            writer.write(
                String.format("%-20d iii (random number seed)\n", iii)
            );
			System.out.print(String.format("%-20d iii (random number seed)\n", iii));
            // Write the length of the sequences.
            writer.write(
                String.format(
                    "%-20d lengthseq (after deleting gaps, etc.)\n",
                    samplePhylogeny.length()
                )
            );
			System.out.print(String.format(
                    "%-20d lengthseq (after deleting gaps, etc.)\n",
                    samplePhylogeny.length()
                ));
            // Write the whichavg value.
            int jwhichavg = masterVariables.getCriterion();
            writer.write(String.format("%-20d jwhichavg\n", jwhichavg));
			System.out.print(String.format("%-20d jwhichavg\n", jwhichavg));
			System.out.println("------------------------------------------------------");
            // Write the probthresh value.
            //writer.write("-0.57700             probthresh\n");
        }
        catch (IOException e) {
            System.out.println("Error writing the input file for the " +
                               "confidence interval program.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the input file for the " +
                                       "confidence interval program.");
                }
            }
        }
    }
	
    protected int [] getDemarcOutput() {
        BufferedReader input = null;
        int [] demarcValues = null;
        try {
            input = new BufferedReader(new FileReader(outputFileName2));
            bestLike = 0.0;
            double nextLike;
            String line = input.readLine();
            StringTokenizer tk;
            demarcValues = new int[3];
            int nextNpop;
            while (line != null) {
                tk = new StringTokenizer(line);
                nextNpop = (int)(new Double(tk.nextToken())).doubleValue();
                tk.nextToken();
                nextLike = -1 * (new Double(tk.nextToken())).doubleValue();
                if (nextLike > bestLike) {
                    bestLike = nextLike;
                    demarcValues[0] = nextNpop;
                }
                line = input.readLine();
            }
            input.close();
            // Read the same file again, only this time to find the confidence interval.
            double probThresh = bestLike / masterVariables.CI_NUMBER;
            input = new BufferedReader(new FileReader(outputFileName2));
            line = input.readLine();
            while (line != null) {
                tk = new StringTokenizer(line);
                nextNpop = (int)(new Double(tk.nextToken())).doubleValue();
                tk.nextToken();
                nextLike = -1 * (new Double(tk.nextToken())).doubleValue();
                // If we have reached the border of the confidence interval, set the value appropriately.
                if (demarcValues[1] == 0) {
                    if (nextLike > probThresh) {
                        demarcValues[1] = nextNpop;
                    }
                }
                // Otherwise check if we are out of the confidence interval yet.
                else if (nextLike < probThresh) {
                    demarcValues[2] = nextNpop - 1;
                    break;
                }
                line = input.readLine();
                if (line == null && ! (nextLike < probThresh))
                    demarcValues[2] = nextNpop;
            }
            input.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return demarcValues;
    }

    /**
     *  Private method to read the output file from the demarcation confidence
     *  interval program.
     *
     *  @param outputFile The file to read from.
     */
    private void readOutputFile(File outputFile) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(outputFile));
            String nextLine = reader.readLine();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer(nextLine);
                // The first line contains the lower value of the confidence
                // interval for demarcation and the likelihood for that
                // value.
                st.nextToken(); // "lower".
                st.nextToken(); // "bound".
                st.nextToken(); // "npop".
                result = new Integer(st.nextToken()).intValue();
                st.nextToken(); // "likelihood".
                likelihood = new Double(st.nextToken()).doubleValue();
                nextLine = reader.readLine();
            }
        }
        catch (IOException e) {
            System.out.println("Error reading the output file from the " +
                               "confidence interval program.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the output file from " +
                                       "the confidence interval program.");
                }
            }
        }
    }


    private String inputFileName;
    private String outputFileName;
    private String inputFileName2;
    private String outputFileName2;
	private double bestLike;
    private MasterVariables masterVariables;
    private Phylogeny phylogeny;
    private Phylogeny samplePhylogeny;
    private Binning binning;
    //private Hillclimb hillclimb;
	private ParameterSet hillclimb;

    private int result = 0;
    private double likelihood = 0.0;

    private int runNumber;

    private boolean hasRun;

}

